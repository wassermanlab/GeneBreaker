import React from 'react';
import InputComp from './inputComp'
import SelectComp from './selectComp'
import SelectList from './selectListComp'
import { host } from '../host'

function GeneralInfo(props) {
  if (props.page !== 1) {
    return null;
  }

  return (
    <React.Fragment>
      <h1>General Info</h1>
      {/* genome */}
      <SelectComp
        title={"Genome Assembly"}
        name={"genome"}
        value={props.genome}
        onChange={props.onChange}
        options={[{ value: "grch38", text: "GRCh38" }, { value: "grch37", text: "GRCh37" }]}/>
      {/* sex */}
      <SelectComp
        title={"Proband Sex"}
        name={"sex"}
        value={props.sex}
        onChange={props.onChange}
        options={[{ value: "XX", text: "XX" }, { value: "XY", text: "XY" }]}/>
      {/* gene symbol */}
      <InputComp
        title={"Gene symbol"}
        name={"gene_name"}
        value={props.gene_name}
        onChange={props.onChange}/>
      {/* transcript list */}
      <SelectList title={"Transcripts"} url={host + 'get_transcripts/' + props.genome + '/' + props.gene_name}
        populatingText={`genome assembly '${props.genome}' and gene symbol '${props.gene_name}'`}
        name={"gene_uid"} value={props.gene_uid} type={"transcript"} onChange={props.onChange}/>
      
    </React.Fragment >
  )
}


export default GeneralInfo;