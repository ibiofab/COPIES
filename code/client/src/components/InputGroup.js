import React from 'react'
import DropdownInput from './DropdownInput'
import TextInput from './TextInput'

const InputGroup = ({field}) => {
    const dropdowns = ["PAM Orientation","Distance Type", "On-target score"];
    const textinputs = ["PAM", "Backbone Sequence", "Restriction enzyme to avoid"];

    if (dropdowns.includes(field)) {
        return (
            <DropdownInput field={field} />
          );
    } else if (textinputs.includes(field)) {
        return (
            <TextInput field={field} />
        );
    } else {
        console.log(field);
        return (
            <div>{field} not functional</div>
        );
    }
  
}

export default InputGroup