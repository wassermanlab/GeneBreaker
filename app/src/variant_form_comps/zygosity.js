import React from 'react';

function Zygosity(props) {
    // XY sex and X or Y gene 
    if (props.sex === "XY" && (props.chrom === "chrY" || props.chrom === "chrX")) {
      return <option value="hemizygous">Hemizygous</option>
    } else if (props.var === 2) { // variant 2
      return <option value="heterozygous">Heterozygous</option>
    } else { // variant 1
      return (
        <React.Fragment>
          <option value="homozygous">Homozygous</option>
          <option value="heterozygous">Heterozygous</option>
        </React.Fragment>)
    }
  }
  

export default Zygosity;