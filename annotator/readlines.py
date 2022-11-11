import os

sum = 0
for file in os.listdir():
    if file.endswith(".py"):
        with open(file, "r") as fh:
            lines = fh.readlines()
            sum += len(lines)

print(sum)
