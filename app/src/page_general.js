import React from 'react';
import {Field, ErrorMessage } from 'formik';

const PageGeneral = (props) => (
  <React.Fragment>
    {/* gene_uid: "",
        genome: "hg38",
        sex: "" */}
    <div>
      <label>Genome</label>
      <Field name="genome" component="select" validate={props.required}>
        <option value="hg38">hg38</option>
        <option value="hg19">hg19</option>
      </Field>
    </div>

    <div>
      <label>Proband Sex</label>
      <Field name="sex" component="select" validate={props.required}>
        <option value="">Select</option>
        <option value="XX">XX</option>
        <option value="XY">XY</option>
      </Field>
      <ErrorMessage
        name="sex"
        component="div"
        className="field-error"
      />
    </div>

    <div>
      <label>Gene UID</label>
      <Field
        name="gene_uid"
        component="input"
        type="text"
        validate={props.required}
      />
      <ErrorMessage
        name="gene_uid"
        component="div"
        className="field-error"
      />
    </div>
  </React.Fragment>

);

export default PageGeneral;
