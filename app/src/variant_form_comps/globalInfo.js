import React from 'react';
import SelectComp from './selectComp'
import SelectList from './selectList'
import host from '../config'

function GlobalInfo(props) {

  if (props.page !== 1) {
    return null;
  }

  return (
    <React.Fragment>
      {/* genome */}
      <SelectComp
        title={"Genome"}
        name={"genome"}
        value={props.genome}
        onChange={props.onChange}
        options={[{ value: "hg38", text: "hg38" }, { value: "hg19", text: "hg19" }]}
      />
      {/* sex */}
      <SelectComp
        title={"Sex"}
        name={"sex"}
        value={props.sex}
        onChange={props.onChange}
        options={[{ value: "XX", text: "XX" }, { value: "XY", text: "XY" }]}
      />
      {/* gene name/symbol */}
      <label>Gene symbol</label>
      <input type="text" className="form-control" name="gene_name" value={props.gene_name} onChange={props.onChange} />
      <SelectList title={"Transcripts"} url={host + 'get_transcripts/' + props.genome + '/' + props.gene_name}
        name={"gene_uid"} value={props.gene_uid} type={"transcript"} onChange={props.onChange}
      />
    </React.Fragment >
  )
}


export default GlobalInfo;