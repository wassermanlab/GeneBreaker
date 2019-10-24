import React, { useState, useEffect } from 'react';
import SelectList from './selectList'
import host from '../config'

function STR(props) {
  const [url, setUrl] = useState(host + 'get_str/' + props.genome + '/' + props.gene_uid + '/' + props.region);

  useEffect(() => {

    function onRegionChange() {
      if (props.region === "CUSTOM") {
        const region = (props.chrom + ":" + props.customStart + "-" + props.customEnd)
        setUrl(host + 'get_str/' + props.genome + '/' + props.gene_uid + '/' + region)
      } else {
        setUrl(host + 'get_str/' + props.genome + '/' + props.gene_uid + '/' + props.region)
      }
    }

    onRegionChange()
  }, [props.chrom, props.customStart, props.customEnd, props.region, props.gene_uid, props.genome]);

  if (props.type !== "str" || props.page - 1 !== props.var) {
    return null;
  }

  return (
    <React.Fragment>
      <SelectList
        title={"ClinGen Variants"}
        url={url}
        name={"str_id_" + props.var}
        value={props.str_id}
        type={"str"} onChange={props.onChange}
      />
      <small className="form-text text-muted">
        Known pathogenic short tandem repeats have a non-zero pathogenicity representing a known repeat difference length.
          </small>
      <label>Repeat difference length</label>
      <div className="input-group" >
        <input type="text" className="form-control" value={props.length} name={"length_" + props.var} onChange={props.onChange} />
      </div>
      <small className="form-text text-muted">
        Positive integers represent the number of additional repeats, negative integers represent the number of retractions.
          </small>
    </React.Fragment>
  )
}

export default STR;