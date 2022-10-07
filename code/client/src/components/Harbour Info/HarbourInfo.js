import React from 'react'
import InputGroup from '../DropdownInput'

const HarbourInfo = () => {
    const fields = [
        "Distance from surrounding genes (bp)",
        "Distance to measure gene density (bp)",
        "Distal end of chromosome to avoid (bp)",
        "Reference organism/s for BLAST"
    ]
  return (
    <div>
        HarbourInfo
        {
            fields.map(field => {
                return (
                    <InputGroup field={field} />
                )
            })
        }

    </div>
  )
}

export default HarbourInfo