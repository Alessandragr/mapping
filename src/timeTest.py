import time
import subprocess
import os  # Import os module to use os.path functions

# Python function to measure
def fileVerifier(filePath):
    """
    Check if the file is a valid SAM file.
    
    :param filePath: path of the file to check
    :return: True if the file is valid, else False
    """
    # Check if it's a regular file
    if not os.path.isfile(filePath):
        print(f"Error: {filePath} is not valid.")
        return False

    # Check if the file has a .sam extension
    if not filePath.endswith('.sam'):
        print(f"Error: {filePath} doesn't have a '.sam' extension.")
        return False

    # Check if the file is empty
    if os.path.getsize(filePath) == 0:
        print(f"Error: file '{filePath}' is empty.")
        return False

    print(f"File '{filePath}' is valid and not empty.")

    # Validate the number of columns in the first 3 non-header lines
    with open(filePath, 'r') as file:
        lineCount = 0
        for line in file:
            # Ignore headers
            if line.startswith('@'):
                continue

            # Count the number of columns separated by tabulation
            numColumns = len(line.strip().split('\t'))

            if numColumns < 11:
                print(f"Error: line '{line.strip()}' only has {numColumns} columns.")
                return False

            lineCount += 1
            
    print(f"File '{filePath}' has the expected number of columns.")
    print(f"The file: '{filePath}' can be used for the next step.")
    return True

# Path to the log file
log_file_path = "/home/najat/mapping/src/timing_results.log"

# Function to log output to a file
def log_output(message):
    with open(log_file_path, 'a') as log_file:
        log_file.write(message + "\n")

# Measure the running time of the Python function
print("running python function")
start_time = time.perf_counter()  # Start time
fileVerifier(filePath="/home/najat/mapping/src/mapping.sam")  # Run the Python function
end_time = time.perf_counter()  # End time
python_time = end_time - start_time  # Calculate the time taken for the Python function
log_output(f"Python function took {python_time:.4f} seconds.")

print("running bash function")
# Bash command to measure (make sure the Bash script logic is inside a file or a command string)
filename = os.path.basename("/home/najat/mapping/src/mapping.sam")  # Extract the filename
bash_command = ["/home/najat/mapping/src/standardFileVerifier.sh", filename]

# Measure the running time of the Bash command
start_time = time.perf_counter()  # Start time
subprocess.run(bash_command, shell=True)  # Run the Bash command
end_time = time.perf_counter()  # End time
bash_time = end_time - start_time  # Calculate the time taken for the Bash command
log_output(f"Bash command took {bash_time:.4f} seconds.")

print(f"Timing results saved to {log_file_path}")
