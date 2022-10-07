import React from 'react'
import InputGroup from '../DropdownInput'

const HomologyArmInfo = () => {
    const fields = [
        "Length of Homology arm (bp)",
        "Restriction enzyme to avoid",
        "PolyG to avoid",
        "PolyT to avoid"
    ]
  return (
    <div>
        HomologyArmInfo
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

export default HomologyArmInfo