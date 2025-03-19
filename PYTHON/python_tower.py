import socket
import json
import os
import subprocess
import time

RECOMPILE = True
AUTO_RUN = True

def send_and_receive_data(host, port, data):
    # Serialize the data to JSON format
    json_data = json.dumps(data)

    print(json_data)

    # Create a socket object
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    # Connect to the server
    s.connect((host, port))

    # Send the JSON data
    s.sendall(json_data.encode('utf-8'))

    # Wait for the response from the server
    response = s.recv(1024).decode('utf-8')

    # Close the connection
    s.close()

    # Print the response
    print(f"Received response: {response}")

    return response

def spin(molecule_json):
    # Define the host and the port
    host = 'localhost'
    port = 12345

    # Compile (dev only)
    if RECOMPILE:
        os.system("g++ -o ./CPP/compiled_cpp/cpp_tower ./CPP/cpp_tower.cpp -lGL -lGLU -lglut")

    if AUTO_RUN:
        # Spin up the shit
        process = subprocess.Popen(["./CPP/compiled_cpp/cpp_tower"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        while True:
            # Read the next line of output
            output = process.stdout.readline().decode()

            # If the output is 'alive\n', break the loop
            if output == 'alive\n':
                break

            # If there is no more output and the process has finished, break the loop
            if output == '' and process.poll() is not None:
                break

            # Sleep for a short time to prevent this loop from using too much CPU
            time.sleep(0.1)

    reponse = send_and_receive_data(host, port, molecule_json)

    r = json.loads(reponse)

    return r