import React, { useState, useEffect } from 'react';
function Str(props) {
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
      const response = await fetch('http://127.0.0.1:5001/get_str/' + props.genome + '/' + props.gene_uid + '/' + props.region);
      const json = await response.json();
      setStr(json);
    }

    fetchUrl();
  }, [props.region, props.gene_uid, props.genome]);

  if (props.type !== "str") {
    return null;
  }

  return (
    <React.Fragment>
      <div className="form-group">
        <label>Short tandem repeat</label>
        <select className="form-control"
          name={"var" + props.var + "_str_id"}
          value={props.str_id}
          onChange={props.handleInputChange}
          size="5">
          <option key="0" value=""></option>
          {str_list.map((item, index) => (
            <option key={index + item.qualifiers.uid} value={item.qualifiers.uid}>
              Start: {item.start + 1} &emsp;&emsp; End: {item.end} &emsp;&emsp; Motif: {item.qualifiers.motif} &emsp;&emsp; Pathogenicity: {item.qualifiers.pathogenicity} </option>
          ))}
        </select>
        <small className="form-text text-muted">
          Known pathogenic short tandem repeats have a non-zero pathogenicity representing a known repeat difference length.
          </small>
      </div>
      <label>Repeat difference length</label>
      <div className="input-group" >
        <input type="text" className="form-control" value={props.length} name={"var" + props.var + "_length"} onChange={props.handleInputChange}/>
      </div>
      <small className="form-text text-muted">
        Positive integers represent the number of additional repeats, negative integers represent the number of retractions.
          </small>
    </React.Fragment>
  )
}


export default Str;