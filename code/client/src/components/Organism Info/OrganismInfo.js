import React from "react";
import Fuse from "fuse.js";
import Select from "react-select";
import { useState } from "react";
const OrganismInfo = () => {
	const list = [
	{value: "Liam", label: "Liam"},
	{value: "Olivia", label: "Olivia"},
	{value: "Noah", label: "Noah"},
	{value: "Emma", label: "Emma"},
	{value: "Oliver", label: "Oliver"},
	{value: "Charlotte", label: "Charlotte"},
	{value: "Elijah", label: "Elijah"},
	{value: "Amelia", label: "Amelia"},
	{value: "James", label: "James"},
	{value: "Ava", label: "Ava"},
	{value: "William", label: "William"},
	{value: "Sophia", label: "Sophia"},
	{value: "Benjamin", label: "Benjamin"},
	{value: "Isabella", label: "Isabella"},
	{value: "Lucas", label: "Lucas"},
	{value: "Mia", label: "Mia"},
	{value: "Henry", label: "Henry"},
	{value: "Evelyn", label: "Evelyn"},
	{value: "Theodore", label: "Theodore"},
	{value: "Harper", label: "Harper"},
	];
	const options = {
		isCaseSensitive: false,
		// includeScore: false,
		// shouldSort: true,
		// includeMatches: false,
		findAllMatches: true,
		// minMatchCharLength: 1,
		// location: 0,
		// threshold: 0.2,
		// distance: 100,
		// useExtendedSearch: false,
		// ignoreLocation: false,
		// ignoreFieldNorm: false,
		// fieldNormWeight: 1,
		// keys: [
		//   "title",
		//   "author.firstName"
		// ]
	  };
	const [state, setState] = useState("");
	const fuse = new Fuse(list, options);

	return <div>OrganismInfo
		<Select options={list} />
		{/* <input placeholder="name" onChange={e=>setState(e.target.value)} />
	  <select>
	  {fuse.search(state).map(element => {
		return (
			<option>{element.item}</option>
		)
	  })}
	  </select> */}

	</div>;
};

export default OrganismInfo;
