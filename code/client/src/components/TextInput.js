import React from "react";

const TextInput = ({ field }) => {
  let regex;
  switch (field) {
    case "PAM":
      regex = new RegExp(/[ACTGRYSWKMBDHVN]/g);
      break;
    case "Backbone Sequence":
      regex = new RegExp(/[ACTG]/g);
      break;
    case "Restriction enzyme to avoid":
      regex = new RegExp(/[ACTGRYSWKMBDHVN]/g);
      break;
    default:
      regex = new RegExp(/[a-zA-Z]/g);
  }
	return (
		<div>
			<div className="input-group">
				<label htmlFor={field}>{field}:</label>
				<input className="input text-input" id={field} pattern={"/[ACTGRYSWKMBDHVN]/g"}/>
			</div>
		</div>
	);
};

export default TextInput;
