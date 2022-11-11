import os

folders = ["results", "log"]
basefolder = os.path.dirname(os.getcwd())

for folder in folders:
    path = os.path.join(basefolder, folder)
    if not os.path.exists(path):
        os.mkdir(path)
