library(tibble)

wikipedia_g <- graph_from_data_frame(
  tibble::tribble(
    ~from, ~to, ~expected_spc, ~expected_main_path
    "A", "C", 2, 0,
    "B", "C", 2, 0,
    "B", "D", 5, 1,
    "B", "J", 1, 0,
    "C", "E", 2, 0,
    "C", "H", 2, 0,
    "D", "F", 3, 1,
    "D", "I", 2, 0,
    "J", "M", 1, 0,
    "E", "G", 2, 0,
    "F", "H", 1, 0,
    "F", "I", 2, 1,
    "G", "H", 2, 0,
    "I", "L", 2, 0,
    "I", "M", 2, 1,
    "H", "K", 5, 0,
    "M", "N", 3, 1),
  directed = TRUE)


all(E(wikipedia_g)$expected_spc == spc(wikipedia_g))
# TRUE
all(E(wikipedia_g)$expected_main_path == main_search(wikipedia_g))