#include <iostream>
#include <string>
#include "../external_tools/json.hpp"  // Make sure the path to the json header is correct for your setup

using json = nlohmann::json;

namespace MDS1 {

json process(const std::string& jsonData) {
    std::string error_msg = "An error occurred during processing";

    try {
        // Parse the JSON data
        auto received_data = json::parse(jsonData);
        
        // Check if the required keys are present and are of the expected type
        if (!received_data.contains("target") || !received_data["target"].is_object()) {
            throw std::runtime_error("Missing or invalid 'target' in JSON");
        }
        if (!received_data.contains("molecules") || !received_data["molecules"].is_object()) {
            throw std::runtime_error("Missing or invalid 'molecules' in JSON");
        }

        // Retrieve the "target" key containing a dict object
        json target = received_data["target"];

        // Save the target object
        std::string targetString = target.dump();

        // Molecule goes here
        json molecule = received_data["molecules"];

        // Save the molecule object
        std::string moleculeString = molecule.dump();

        // Prepare a response json
        json response = {{"target", targetString}, {"molecules", moleculeString}, {"result", "Success"}};

        // Return the response JSON object
        return response;
    } catch (json::parse_error& e) {
        // Handle JSON parsing errors (e.g., invalid JSON)
        std::cerr << "JSON Parsing Error: " << e.what() << '\n';
        error_msg = "Invalid JSON data 1";
    } catch (json::type_error& e) {
        // Handle JSON type errors (e.g., accessing a key that doesn't exist or is of an unexpected type)
        std::cerr << "JSON Type Error: " << e.what() << '\n';
        error_msg = "Invalid JSON data 2";
    } catch (std::exception& e) {
        // Handle any other std::exception derived errors
        std::cerr << "Error: " << e.what() << '\n';
        error_msg = e.what();
    } catch (...) {
        // Handle any other errors
        std::cerr << "An unknown error occurred\n";
        error_msg = "An unknown error occurred";
    }

    // Return an error response if any exception is caught
    return json{{"error", "An error occurred during processing: " + error_msg}};
}

} // namespace MDS1
