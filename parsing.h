

bool line_not_readable(const std::string &line)
{
    return (line[0] == '#' || line.empty());
}

std::string clear_line(const std::string &initial_line)
{
    std::string line = initial_line;
    line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
                   return std::isspace(static_cast<unsigned char>(c));
               }),
               line.end());
    return line;
}

void create_boundary_forces(const std::string &initial_data, double *Forces, Node *DNodes, long node_index)
{
    std::string data = initial_data;
    int delimiterPos = 0;
    delimiterPos = data.find(",");
    for (int DOF_Node_index = 0; DOF_Node_index < DNodes[node_index].DOF; DOF_Node_index++)
    {
        std::string value = data.substr(0, delimiterPos);
        data = data.substr(delimiterPos + 1);
        Forces[DNodes[node_index].absolute_DOF_addr + DOF_Node_index] = std::atol(value.c_str());
        delimiterPos = data.find(",");
    }
}

void create_boundary_displacement(const std::string &initial_data, int *indx, Node *DNodes, long node_index)
{
    std::string data = initial_data;
    int delimiterPos = 0;
    delimiterPos = data.find(",");
    for (int DOF_Node_index = 0; DOF_Node_index < DNodes[node_index].DOF; DOF_Node_index++)
    {
        std::string value = data.substr(0, delimiterPos);
        data = data.substr(delimiterPos + 1);
        indx[DNodes[node_index].absolute_DOF_addr + DOF_Node_index] = std::stoi(value, nullptr, 10); // indx cest les contraintes ?
        delimiterPos = data.find(",");
    }
}