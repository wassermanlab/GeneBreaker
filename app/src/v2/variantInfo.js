import React, { useState, useEffect } from 'react';
import VariantDetails from './variantDetails'
import SelectComp from './selectComp'
import { host } from '../host'

function VariantInfo(props) {
  const [url, setUrl] = useState("");
  useEffect(() => {

    function onRegionChange() {
      let region
      if (props.region === "CUSTOM") {
        region = (props.chrom + ":" + props.customStart + "-" + props.customEnd)
      } else { region = props.region }
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
  if (props.page === 2) {
    if (props.sex === "XY" && (props.chrom === "chrX" || props.chrom === "chrY")) {
      zygosity_options = [{ value: "", text: "Select" }, { value: "HEMIZYGOUS", text: "hemizygous" }];
    } else {
      zygosity_options.push({ value: "HOMOZYGOUS", text: "homozygous" })
    }
  }
  if (props.page < 2 || props.page > 3) {
    return null;
  }
  else if (props.page === 3 && props.zygosity_1 !== "HETEROZYGOUS") {
    return (
      <React.Fragment>
        <p>Variant 2 cannot be made with a homozygous or hemizygous variant 1.</p>
      </React.Fragment>
    )
  }
  return (
    <React.Fragment>
      {/*  */}
      {/* Region */}
      <SelectComp
        title={"Region"}
        name={"region_" + (props.page - 1)}
        value={props.region}
        onChange={props.onChange}
        options={[{ value: "", text: "Select" }, { value: "CODING", text: "coding" }, { value: "GENIC", text: "genic" },
        { value: "UTR", text: "Untranslate region" }, { value: "INTRONIC", text: "intronic" }, { value: "CUSTOM", text: "custom" }]} />
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
     url={url} type={props.type} var={props.page - 1}
     clinvar_id={props.clinvar_id_1} start={props.start_1} end={props.end_1}
     clingen_id={props.clingen_id_1} length={props.length_1} element={props.element_1}
     snv_type={props.snv_type_1} str_id={props.str_id_1} onChange={props.handleInputChange}
    />
    </React.Fragment >
  )
}


export default VariantInfo;