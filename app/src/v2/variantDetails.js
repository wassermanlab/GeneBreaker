import React from 'react';
import InputComp from './inputComp'
import SelectList from './selectListComp'
import SelectComp from './selectComp'

function VariantInfo(props) {
  switch (props.type) {
    case "":
      return null;
    case "CLINVAR":
      return (<SelectList title={"Clinvar variants"} url={props.url} name={"clinvar_id_" + props.var}
        value={props.clinvar_id} type={"clinvar"} onChange={props.onChange} />)
    case "CLINGEN":
      return (<SelectList title={"ClinGen variants"} url={props.url} name={"clingen_id_" + props.var}
        value={props.clingen_id} type={"clingen"} onChange={props.onChange} />)
    case "CNV":
      return (<React.Fragment>
        <InputComp title={"Start position"} name={"start_" + props.var} value={props.start} onChange={props.onChange} />
        <InputComp title={"End position"} name={"end_" + props.var} value={props.end} onChange={props.onChange} />
        <SelectComp
          title={"Copy change"}
          name={"length_" + props.var}
          value={props.length}
          onChange={props.onChange}
          options={[{ value: "", text: "Select" }, { value: -1, text: "-1" }, { value: 1, text: "1" }]} />
      </React.Fragment>)
    case "MEI":
      return (<React.Fragment>
        <InputComp title={"Start position"} name={"start_" + props.var} value={props.start} onChange={props.onChange} />
        <SelectComp
          title={"Element type"}
          name={"element_" + props.var}
          value={props.element}
          onChange={props.onChange}
          options={[{ value: "", text: "Select" }, { value: "ALU", text: "ALU" },
          { value: "LINE", text: "LINE" }, { value: "SINE", text: "SINE" }]} />
      </React.Fragment>)
    case "INDEL":
      return (<React.Fragment>
        <InputComp title={"Start position"} name={"start_" + props.var} value={props.start} onChange={props.onChange} />
        <InputComp title={"Length"} name={"length_" + props.var} value={props.length} onChange={props.onChange} />
      </React.Fragment>)
    case "SNV":
      return (<React.Fragment>
        <InputComp title={"Start position"} name={"start_" + props.var} value={props.start} onChange={props.onChange} />
        <SelectComp
          title={"SNV type"}
          name={"snv_type_" + props.var}
          value={props.snv_type}
          onChange={props.onChange}
          options={[{ value: "", text: "Select" }, { value: "STOPLOSS", text: "stoploss" }, { value: "MISSENSE", text: "missense" },
          { value: "NONSENSE", text: "nonsense" }, { value: "SYNONYMOUS", text: "synonymous" }, { value: "A", text: "A" },
          { value: "T", text: "T" }, { value: "C", text: "C" }, { value: "G", text: "G" }, { value: "ANY", text: "Any" }]} />
      </React.Fragment>)
    case "STR":
      return (<React.Fragment>
        <SelectList title={"Short tandem repeats"} url={props.url} name={"str_id_" + props.var}
          value={props.str_id} type={"str"} onChange={props.onChange} />
        <InputComp title={"Repeat length"} name={"length_" + props.var} value={props.length} onChange={props.onChange} />
      </React.Fragment>)
    default:
      return null;
  }

}


export default VariantInfo;