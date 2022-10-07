import React from "react";
import InputGroup from "../InputGroup";

const GuideRNAInfo = () => {
	const fields = [
		"PAM",
		"PAM Orientation",
		"GC content",
		"Guide Length",
		"Seed Length",
		"PolyG to avoid",
		"Distance Type",
		"On-target score",
		"PolyT to avoid",
		"Edit Distance",
		"Backbone Sequence",
		"Restriction enzyme to avoid",
	];
	return (
		<div>
			GuideRNAInfo
			{fields.map((field) => {
				return (
					<>
						<InputGroup field={field} />
					</>
				);
			})}
		</div>
	);
};

export default GuideRNAInfo;
