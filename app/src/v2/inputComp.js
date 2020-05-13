import React from 'react';

function InputComp(props) {
  return (
    <React.Fragment>
      <div className="form-group row">
        <label className="col-sm-2 col-form-label">{props.title}</label>
        <div className="col-sm-10">
          <input type="text" className="form-control" name={props.name} value={props.value} onChange={props.onChange} />
          <small class="form-text text-muted">
          {props.help}
        </small>
        </div>
      </div>
    </React.Fragment>
  )
}
export default InputComp;
