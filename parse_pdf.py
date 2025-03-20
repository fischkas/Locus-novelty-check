import fitz  # PyMuPDF
import pandas as pd

# File paths
keywords_file = "novelty_lookup_search_file.txt"
pdf_file = "scientific_articles_daymonthyear.txt"
output_file = "novelty_lookup_genes.results"

# Step 1: Load keywords
with open(keywords_file, "r", encoding="utf-8") as f:
    keywords = [line.strip() for line in f.readlines()]

# Step 2: Extract text from the PDF and track page numbers
doc = fitz.open(pdf_file)
keyword_results = []


for page_num, page in enumerate(doc):
    page_text = page.get_text("text").lower()
    for keyword in keywords:
        count = page_text.count(keyword.lower())
        if count > 0:
            keyword_results.append((keyword, count, page_num + 1))  # Store keyword, count, and page number

# Step 3: Convert to DataFrame and save results
df = pd.DataFrame(keyword_results, columns=["Keyword", "Occurrences", "Page Number"])
df.to_csv(output_file, index=False)

print(f"Keyword search completed. Results saved to {output_file}.")
