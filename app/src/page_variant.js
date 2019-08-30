import React from 'react';
import {Field, ErrorMessage } from 'formik';

const PageVariant = (props) => (
  <React.Fragment>
    {console.log(props)}
    <div>
      <label>Zygosity</label>
      <Field name={props.var+".zygosity"} component="select" validate={props.required}>
        <option value="">Select</option>
        <option value="hemizygous">Hemizygous</option>
        <option value="heterozygous">Heterozygous</option>
        <option value="homorozygous">Homorozygous</option>
      </Field>
      <ErrorMessage
        name={props.var+".zygosity"}
        component="div"
        className="field-error"
      />
    </div>
    <div>
      <label>Type</label>
      <Field name={props.var+".type"} component="select" validate={props.required}>
        <option value="">Select</option>
        <option value="clinvar">ClinVar</option>
        <option value="clingen_cnv">ClinGen Copy Number Variant</option>
        <option value="cnv">Copy Number Variant</option>
        <option value="indel">Indel</option>
        <option value="mei">Mobile Element Insertion</option>
        <option value="snv">Single Nucleotide Variant</option>
        <option value="str">Short Tandem Repeat</option>
      </Field>
      <ErrorMessage
        name={props.var+".type"}
        component="div"
        className="field-error"
      />
    </div>
    <div>
      <label>Region</label>
      <Field name={props.var+".region"} component="select" validate={props.required}>
        <option value="">Select</option>
        <option value="coding">Coding</option>
        <option value="utr">Untranslated Region</option>
        <option value="intronic">Intronic</option>
        <option value="genic">Genic</option>
        <option value="custom">Custom</option>
      </Field>
      <ErrorMessage
        name={props.var+".region"}
        component="div"
        className="field-error"
      />
    </div>
  </React.Fragment>

);

export default PageVariant;
