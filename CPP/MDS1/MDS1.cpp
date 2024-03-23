#include <iostream>
#include <string>
#include "../external_tools/json.hpp"  // Make sure the path to the json header is correct for your setup

using json = nlohmann::json;

namespace MDS1 {

// Function to process the input JSON and return a response JSON
json process(const std::string& jsonData) {
    // Parse the JSON data
    auto received_data = json::parse(jsonData);
    

    // Target goes here

    // Molecule goes here



    // For demonstration, let's just prepare a response json
    json response = {{"message", "Data (BLOOOOOOOOOO) received successfully"}};

    // Return the response JSON object
    return response;
}

} // namespace MDS
