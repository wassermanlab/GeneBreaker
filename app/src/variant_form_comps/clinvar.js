import React, { useState, useEffect } from 'react';
import SelectList from './selectList'

function Clinvar(props) {
  const [url, setUrl] = useState('http://127.0.0.1:5001/get_clinvar/' + props.genome + '/' + props.gene_uid + '/' + props.region);

  useEffect(() => {

    function onRegionChange() {
      if (props.region === "CUSTOM") {
        const region = (props.chrom + ":" + props.customStart + "-" + props.customEnd)
        setUrl('http://127.0.0.1:5001/get_clinvar/' + props.genome + '/' + props.gene_uid + '/' + region)
      } else {
        setUrl('http://127.0.0.1:5001/get_clinvar/' + props.genome + '/' + props.gene_uid + '/' + props.region)
      }
    }

    onRegionChange()
  }, [props.chrom, props.customStart, props.customEnd, props.region, props.gene_uid, props.genome]);

  if (props.type !== "clinvar" || props.page - 1 !== props.var) {
    return null;
  }

  return (
    <React.Fragment>
      <SelectList
        title={"ClinVar Variants"}
        url={url}
        name={"clinvar_id_" + props.var}
        value={props.clinvar_id}
        type={"clinvar"} onChange={props.onChange}
      />
    </React.Fragment>
  )
}

export default Clinvar;