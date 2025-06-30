require "./library/graph"
require "./library/utilities"

# Check expected CLI arguments are present
if ARGV.size != 3
  puts "Error: Unexpected number of arguments"
  puts "Usage: ./network [input_file] [output_file] [k]"
  exit 1
end

# Parse CLI arguments
input_file = ARGV[0]
output_file = ARGV[1]
k = ARGV[2].to_i

# Construct network and write contents to the specified file
graph, error = read_data(input_file)
unless error == ""
  puts error
  exit 1
end

network = graph.construct_network(k)

save_graph(output_file, network)
