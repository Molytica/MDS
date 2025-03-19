#include <iostream>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <string>
#include "external_tools/json.hpp"  // Make sure this path is correct based on your setup
#include "MDS1/MDS1.cpp"
#include "VIZ/VIZ.cpp"

using json = nlohmann::json;

// Function modified to take a JSON string and return a JSON object
json processAndRespond(const std::string& jsonData) {
    
    // Process the json data input and create the json response
    VIZ::visualize(jsonData);
    
    return MDS1::process(jsonData);
}


int main() {
    // Print "alive" to the terminal output
    std::cout << "alive" << std::endl;

    int server_fd, new_socket;
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);
    char buffer[1024] = {0};

    // Creating socket file descriptor
    if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
        perror("socket failed");
        exit(EXIT_FAILURE);
    }

    // Forcefully attaching socket to the port 12345
    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt))) {
        perror("setsockopt");
        exit(EXIT_FAILURE);
    }
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(12345);

    // Forcefully attaching socket to the port 12345
    if (bind(server_fd, (struct sockaddr *)&address, sizeof(address)) < 0) {
        perror("bind failed");
        exit(EXIT_FAILURE);
    }
    if (listen(server_fd, 3) < 0) {
        perror("listen");
        exit(EXIT_FAILURE);
    }
    if ((new_socket = accept(server_fd, (struct sockaddr *)&address, (socklen_t*)&addrlen)) < 0) {
        perror("accept");
        exit(EXIT_FAILURE);
    }

    // Read the data sent by the client
    read(new_socket, buffer, 1024);

    // Process the received data and get the response
    json response = processAndRespond(std::string(buffer));

    // Send the response back to the client
    send(new_socket, response.dump().c_str(), response.dump().length(), 0);

    // Close the sockets
    close(new_socket);
    close(server_fd);

    return 0;
}
