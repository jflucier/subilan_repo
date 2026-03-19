# 1. Install/Load required libraries
if (!require("arrow")) install.packages("arrow")
if (!require("reshape2")) install.packages("reshape2")

library(arrow)
library(reshape2)

# 2. Load the parquet report
# Replace 'report.parquet' with your actual file path
report <- read_parquet("report.parquet")

# 3. Apply standard DIA-NN filters for high-confidence hits
# Note: These are typical thresholds for a 1% FDR pipeline
filtered_data <- report[
  report$Q.Value <= 0.01 &
  report$PG.Q.Value <= 0.01 &
  !is.na(report$PG.MaxLFQ),
]

# 4. Extract only the necessary columns to save memory
# We need 'Run' (samples), 'Protein.Group' (rows), and 'PG.MaxLFQ' (values)
simplified_df <- filtered_data[, c("Run", "Protein.Group", "PG.MaxLFQ")]

# 5. Pivot from "Tall" format to "Wide" Matrix
# This uses MaxLFQ intensity values
pg_matrix <- dcast(
  simplified_df,
  Protein.Group ~ Run,
  value.var = "PG.MaxLFQ",
  fun.aggregate = max # Handle potential duplicates
)

# 6. Save the result
write.table(pg_matrix, "custom_pg_matrix.tsv", sep="\t", row.names = FALSE, quote = FALSE)

# Check the first few rows
head(pg_matrix)
