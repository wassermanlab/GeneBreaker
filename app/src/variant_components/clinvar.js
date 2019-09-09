import React, { useState, useEffect } from 'react';
function Clinvar(props) {
  ///get_str/<genome>/<transcript_uid>/<region>
  const [str_list, setStr] = useState([]);

  function check_region() {
    const re = new RegExp('^chr(1[0-9]|Y|Z|2[0-2]|[1-9])\:\d+\-\d+$/')
    if (props.region in ["CODING", "GENIC", "UTR", "INTRONIC"]) {
      return true;
    } else {
      if (re.test(props.region)) {
        const start = props.regions.split(":")[1].split("-")[0]
        const end = props.regions.split(":")[1].split("-")[0]
        if (parseInt(end) > parseInt(start)) {
          return true;
        }
      }
    }
    return false;
  }

  async function fetchUrl() {
    if (!check_region) {
      return null;
    }
    const response = await fetch('http://127.0.0.1:5001/get_clinvar/' + props.genome + '/' + props.gene_uid + '/' + props.region);
    const json = await response.json();
    setStr(json);
  }

  useEffect(() => {
    fetchUrl();
  }, [props.region]);

  if (props.type !== "clinvar") {
    return null;
  }

  return (
    <React.Fragment>
      <div className="form-group">
        <label>ClinVar variants</label>
        <select className="form-control"
          name={"var" + props.var + "_clinvar_id"}
          value={props.clinvar_id}
          onChange={props.handleInputChange}
          size="5">
          <option key="0" value=""></option>
          {str_list.map((item, index) => (
            <option key={item.qualifiers.uid} value={item.qualifiers.uid}>
              Start: {item.start + 1} &emsp;&emsp; Ref: {item.qualifiers.ref} &emsp;&emsp; Alt: {item.qualifiers.alt} &emsp;&emsp; Clinical Significance: {item.qualifiers.CLNSIG} </option>
          ))}
        </select>
      </div>
    </React.Fragment>
  )
}


export default Clinvar;