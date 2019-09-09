import React, { useState, useEffect } from 'react';
function Clingen(props) {
  ///get_str/<genome>/<transcript_uid>/<region>
  const [str_list, setStr] = useState([]);



  useEffect(() => {
    function check_region() {
      // eslint-disable-next-line
      const re = new RegExp('^chr(1[0-9]|Y|Z|2[0-2]|[1-9])\:\d+\-\d+$/')
      if (props.region in ["CODING", "GENIC", "UTR", "INTRONIC"]) {
        return true;
      } else {
        if (re.test(props.region)) {
          const start = props.region.split(":")[1].split("-")[0]
          const end = props.region.split(":")[1].split("-")[0]
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
      const response = await fetch('http://127.0.0.1:5001/get_clingen/' + props.genome + '/' + props.gene_uid + '/' + props.region);
      const json = await response.json();
      setStr(json);
    }

    fetchUrl();
  }, [props.region, props.gene_uid, props.genome]);

  if (props.type !== "clingen") {
    return null;
  }

  return (
    <React.Fragment>
      <div className="form-group">
        <label>GlinGen variants</label>
        <select className="form-control"
          name={"var" + props.var + "_clingen_uid"}
          value={props.clingen_uid}
          onChange={props.handleInputChange}
          size="5">
          <option key="0" value=""></option>
          {str_list.map((item, index) => (
            <option key={index + item.qualifiers.uid} value={item.qualifiers.uid} start={item.start + 1} end={item.end} copy_change={item.qualifiers.copy_number_change}>
              Start: {item.start + 1} &emsp;&emsp; End: {item.end} &emsp;&emsp; Clinical assertion: {item.qualifiers.clinical_assertion} &emsp;&emsp; Copy change: {item.qualifiers.copy_number_change} </option>
          ))}
        </select>
      </div>
    </React.Fragment>
  )
}

export default Clingen;