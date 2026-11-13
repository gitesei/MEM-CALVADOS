import json
import numpy as np

with open("Grb2.json") as f:
    data = json.load(f)

pae = data[0]["predicted_aligned_error"]

cropped = [row[:154] for row in pae[:154]]

new_data = [{"predicted_aligned_error": cropped,
            "max_predicted_aligned_error": data[0]["max_predicted_aligned_error"]}]

with open("Grb2DSH3.json", "w") as f:
    json.dump(new_data, f)
