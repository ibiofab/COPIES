import React from "react";
import Dropdown from "react-bootstrap/Dropdown";

const TextInput = ({ field }) => {
  let elements;
  switch (field) {
    case "PAM Orientation":
      elements = ["5prime", "3prime"];
      break;
    case "Distance Type":
      elements = ["hamming", "levenshtein"];
      break;
    case "On-target score":
      elements =["doench", "crospr", "deepguide(Cas9)","deepguide(Cas12a)", "sgRNA_ecoli(Cas9)", "sgRNA_ecoli(eSpCas9"];
      break;
    default:
      elements = ["NOT A DROPDOWN"];
  }
	return (
		<div>
			<div className="input-group">
				<label htmlFor={field}>{field}:</label>
				<select className="input dropdown-input" id={field}>
          {elements.map((element) => {
            return (
              <option>{element}</option>
            );
          })}
        </select>
			</div>
		</div>
	);
};

export default TextInput;
