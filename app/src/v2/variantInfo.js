import React, { useState, useEffect } from 'react';
import VariantDetails from './variantDetails'
import SelectComp from './selectComp'
import { host } from '../host'

function VariantInfo(props) {
  const [url, setUrl] = useState("");
  const [region, setRegion] = useState("");
  useEffect(() => {

    function onRegionChange() {
      let region
      if (props.region === "CUSTOM") {
        region = (props.chrom + ":" + props.customStart + "-" + props.customEnd)
      } else { 
        region = props.region 
      }
      setRegion(props.region)
      switch (props.type) {
        case "CLINVAR":
          setUrl(host + 'get_clinvar/' + props.genome + '/' + props.gene_uid + '/' + region);
          break;
        case "CLINGEN":
          setUrl(host + 'get_clingen/' + props.genome + '/' + props.gene_uid + '/' + region);
          break;
        case "STR":
          setUrl(host + 'get_str/' + props.genome + '/' + props.gene_uid + '/' + region);
          break;
        default:
          break;
      }
    }
    onRegionChange()
  }, [props.type, props.chrom, props.customStart, props.customEnd, props.region, props.gene_uid, props.genome]);

  let zygosity_options = [{ value: "", text: "Select" }, { value: "HETEROZYGOUS", text: "heterozygous" }];
  if (props.var === 1) {
    if (props.sex === "XY" && (props.chrom === "X" || props.chrom === "Y")) {
      zygosity_options = [{ value: "", text: "Select" }, { value: "HEMIZYGOUS", text: "hemizygous" }];
    } else {
      zygosity_options.push({ value: "HOMOZYGOUS", text: "homozygous" })
    }
  }
  if (props.page-1 !== props.var) {
    return null;
  }
  else if (props.var === 2 && (props.zygosity_1 !== "HETEROZYGOUS")) {
    return (
      <React.Fragment>
        <p>Variant 2 cannot be made with a homozygous or hemizygous variant 1.</p>
      </React.Fragment>
    )
  }
  if (props.render === "true") {
  return (
    <React.Fragment>
      <h1>Variant {props.var} Info</h1>
      <SelectComp
        title={"Region"}
        name={"region_" + (props.page - 1)}
        value={props.region}
        onChange={props.onChange}
        options={[{ value: "", text: "Select" }, { value: "CODING", text: "coding" }, { value: "GENIC", text: "genic" },
        { value: "UTR", text: "untranslate region (UTR)" }, { value: "INTRONIC", text: "intronic" }, { value: "CUSTOM", text: "custom" }]} />
      {/* custom regions */}
      {props.region === "CUSTOM" &&
        <div className="form-group row">
          <label className="col-sm-2 col-form-label">Custom region</label>
          <div className="input-group col-sm-10">
            <div className="input-group-prepend">
              <span className="input-group-text">{props.chrom}</span>
            </div>
            <input type="text" name={"customStart_" + (props.page - 1)} value={props.customStart} onChange={props.onChange} className="form-control" />
            <input type="text" name={"customEnd_" + (props.page - 1)} value={props.customEnd} onChange={props.onChange} className="form-control" />
          </div>
        </div>
      }
      {/* Type */}
      <SelectComp
        title={"Type"}
        name={"type_" + (props.page - 1)}
        value={props.type}
        onChange={props.onChange}
        options={[{ value: "", text: "Select" }, { value: "CLINVAR", text: "clinvar" }, { value: "CLINGEN", text: "clingen copy number variant" },
        { value: "CNV", text: "copy number variant" }, { value: "MEI", text: "mobile element insertion" }, { value: "INDEL", text: "indel" },
        { value: "SNV", text: "single nucleotide variant" }, { value: "STR", text: "short tantem repeat" }]} />
      {/* Zygosity */}
      <SelectComp
        title={"Zygosity"}
        name={"zygosity_" + (props.page - 1)}
        value={props.zygosity}
        onChange={props.onChange}
        options={zygosity_options} />
    {/* specif variant details */}
    <VariantDetails  
     url={url} region={region} type={props.type} var={props.var}
     clinvar_id={props.clinvar_id} start={props.start} end={props.end}
     clingen_id={props.clingen_id} length={props.length} element={props.element}
     snv_type={props.snv_type} str_id={props.str_id} onChange={props.onChange}
    />
    
    </React.Fragment >
  )}
  return null;
}


export default VariantInfo;