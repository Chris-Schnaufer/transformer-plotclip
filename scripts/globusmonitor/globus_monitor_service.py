#!/usr/bin/python

import os, shutil, json, time, datetime
from functools import wraps
from flask import Flask, request, Response
from flask.ext import restful
from flask_restful import reqparse, abort, Api, Resource
from globusonline.transfer.api_client import TransferAPIClient, goauth


config = {}
configFile = "config.json"

"""Active task object is of the format:
[{
    "user":                     globus username
    "globus_id":                globus job ID of upload
    "files":        [{          list of files included in task, each with
            "path": "file1",        ...file path, which is updated with path-on-disk once completed
            "md": {}                ...metadata to be associated with that file
        }, {...}, ...],
    "received":                 timestamp when task was sent to monitor API
    "completed":                timestamp when task was completed (including errors and cancelled tasks)
    "status":                   can be "IN PROGRESS", "DONE", "ABORTED", "ERROR"
}, {...}, {...}, ...]"""
activeTasks = {}

app = Flask(__name__)
api = restful.Api(app)


# ----------------------------------------------------------
# SHARED UTILS
# ----------------------------------------------------------
"""Load contents of .json file into a JSON object"""
def loadJsonFile(filename):
    f = open(filename)
    jsonObj = json.load(f)
    f.close()
    return jsonObj

"""Load activeTasks from file into memory"""
def loadActiveTasksFromDisk():
    activePath = config.api.activePath

    # Prefer to load from primary file, try to use backup if primary is missing
    if not os.path.exists(activePath):
        if os.path.exists(activePath+".backup"):
            shutil.copyfile(activePath+".backup", activePath)
        else:
            # Create an empty file if primary+backup don't exist
            f = open(activePath, 'w')
            f.write("{}")
            f.close()

    activeTasks = loadJsonFile(activePath)

"""Write activeTasks from memory into file"""
def writeActiveTasksToDisk():
    # Write current file to backup location before writing current file
    activePath = config.api.activePath
    shutil.move(activePath, activePath+".backup")
    f = open(activePath, 'w')
    f.write(json.dumps(activeTasks))
    f.close()

"""Write a completed task onto disk in appropriate folder hierarchy"""
def writeCompletedTaskToDisk(task):
    completedPath = config.api.completedPath
    taskID = task.globus_id

    # e.g. TaskID "eaca1f1a-d400-11e5-975b-22000b9da45e"
    #   = <completedPath>/ea/ca/1f/1a/eaca1f1a-d400-11e5-975b-22000b9da45e.json

    # Create root directory if necessary
    if not os.path.exists(completedPath):
        os.mkdir(completedPath)
    # Create nested hierarchy folders if needed, to hopefully avoid a long flat list
    treeLv1 = os.path.join(completedPath, taskID[:2])
    treeLv2 = os.path.join(treeLv1, taskID[2:4])
    treeLv3 = os.path.join(treeLv2, taskID[4:6])
    treeLv4 = os.path.join(treeLv3, taskID[6:8])
    for dir in [treeLv1, treeLv2, treeLv3, treeLv4]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    # Write to json file with task ID as filename
    dest = os.path.join(treeLv4, taskID+".json")
    f = open(dest, 'w')
    f.write(json.dumps(task))
    f.close()

"""Return full path name to completed logfile for a given task id if it exists, otherwise None"""
def getCompletedTaskLogPath(taskID):
    completedPath = config.api.completedPath

    treeLv1 = os.path.join(completedPath, taskID[:2])
    treeLv2 = os.path.join(treeLv1, taskID[2:4])
    treeLv3 = os.path.join(treeLv2, taskID[4:6])
    treeLv4 = os.path.join(treeLv3, taskID[6:8])
    fullPath = os.path.join(treeLv4, taskID+".json")

    if os.path.isfile(fullPath):
        return fullPath
    else:
        return None

# ----------------------------------------------------------
# API COMPONENTS
# ----------------------------------------------------------
"""Authentication components for API (http://flask.pocoo.org/snippets/8/)"""
def check_auth(username, password):
    """Called to check whether username/password is valid"""
    if username in config.globus.validUsers.keys:
        return password == config.globus.validUsers[username].password
    else:
        return False

def authenticate():
    """Send 401 response that enables basic auth"""
    return Response("Could not authenticate. Please provide valid credentials",
                    401, {"WWW-Authenticate": 'Basic realm="Login Required"'})

def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth or not check_auth(auth.username, auth.password):
            return authenticate()
        return f(*args, **kwargs)
    return decorated

"""Post new globus tasks to be monitored"""
@requires_auth
class GlobusMonitor(restful.Resource):

    """Return list of all active tasks"""
    def get(self):
        return activeTasks, 200

    """Add new Globus task ID from a known user for monitoring"""
    def post(self):
        json_data = request.get_json(force=True)
        for task in json_data:
            # Add to active list if globus username is known, and write to disk
            if task.user in config.globus.validUsers.keys:
                activeTasks[task.globus_id] = {
                    "user": task.user,
                    "globus_id": task.globus_id,
                    "files": task.files,
                    "received": datetime.datetime.now(),
                    "completed": None,
                    "status": "IN PROGRESS"
                }
                writeActiveTasksToDisk()

        return 201

"""Get status of a particular task by globus id"""
@requires_auth
class GlobusTask(restful.Resource):

    """Check if the Globus task ID is finished, in progress, or an error has occurred"""
    def get(self, globusID):
        if globusID not in activeTasks.keys:
            # If not in active list, check for record of completed task
            logPath = getCompletedTaskLogPath(globusID)
            if logPath:
                return loadJsonFile(logPath), 200
            else:
                return "Globus ID not found in active or completed tasks", 404
        else:
            return activeTasks[globusID], 200

    """Remove task from active tasks"""
    def delete(self, globusID):
        if globusID in activeTasks.keys:
            # TODO: Should this allow deletion within Globus as well? For now, just deletes from monitoring
            task = activeTasks[globusID]

            # Write task as completed with an aborted status
            task.status = "DELETED"
            task.completed = datetime.datetime.now()
            task.path = ""

            writeCompletedTaskToDisk(task)
            del activeTasks[task.globus_id]
            writeActiveTasksToDisk()
            return 204

# Add a new Globus id that should be monitored
api.add_resource(GlobusMonitor, '/tasks')
# Check to see if Globus id is finished
api.add_resource(GlobusTask, '/tasks/<string:globusID>')

# ----------------------------------------------------------
# SERVICE COMPONENTS
# ----------------------------------------------------------
"""Use globus goauth tool to get access tokens for valid accounts"""
def generateAuthTokens():
    for validUser in config.globus.validUsers.keys:
        config.globus.validUsers[validUser].authToken = goauth.get_access_token(
                username=validUser,
                password=config.globus.validUsers[validUser].password
            ).token

"""Query Globus API to get current transfer status of a given task"""
def checkGlobusStatus(task):
    authToken = config.globus.validUsers[task.user].authToken

    api = TransferAPIClient(username=task.user, goauth=authToken)
    status_code, status_message, task_data = api.task(task.globus_id)
    if status_code == 200:
        return task_data.status
    else:
        return "UNKNOWN ("+status_code+": "+status_message+")"

"""Continually check globus API for task updates"""
def globusMonitorLoop():
    refreshTimer = 0
    while True:
        time.sleep(1)

        # For tasks whose status is still in-progress, check Globus for transfer status updates
        for task in activeTasks:
            globusStatus = checkGlobusStatus(task)

            if globusStatus in ["SUCCEEDED", "FAILED"]:
                # Update task parameters
                task.status = globusStatus
                task.completed = datetime.datetime.now()
                task.path = getCompletedTaskLogPath(task)

                # Notify Clowder to process file if transfer successful
                if globusStatus == "SUCCEEDED":
                    notifyClowderOfCompletedTask(task)

                # Write out results file, then delete from active list and write new active file
                writeCompletedTaskToDisk(task)
                del activeTasks[task.globus_id]
                writeActiveTasksToDisk()

        # Refresh auth tokens every 12 hours
        refreshTimer += 1
        if refreshTimer >= 43200:
            generateAuthTokens()
            refreshTimer = 0

"""Send Clowder necessary details to load local file after Globus transfer complete"""
def notifyClowderOfCompletedTask(task):
    pass


if __name__ == '__main__':
    loadJsonFile(configFile)
    loadActiveTasksFromDisk()
    generateAuthTokens()

    # Create thread for API to begin listening - requires valid Globus user/pass
    thread.start_new_thread(app.run, kwargs={
        "host": "0.0.0.0",
        "port": int(config.api.port),
        "debug": True
    })
    print("API now listening on port "+config.api.port)

    # Create thread for service to begin monitoring
    thread.start_new_thread(globusMonitorLoop)
    print("Service now monitoring Globus tasks")
