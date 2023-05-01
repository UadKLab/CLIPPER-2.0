import os

# Create directories for results and logs if they don't exist
base_folder = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
result_folder = f"{base_folder}/results"
log_folder = f"{base_folder}/log"

if not os.path.exists(result_folder):
    os.makedirs(result_folder)

if not os.path.exists(log_folder):
    os.makedirs(log_folder)