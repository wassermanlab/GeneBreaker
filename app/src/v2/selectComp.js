import React from 'react';

function SelectComp(props) {
  return (
    <React.Fragment>
      <div className="form-group row">
        <label className="col-sm-2 col-form-label">{props.title}</label>
        <div className="col-sm-10">
          <select
            className="form-control"
            name={props.name}
            value={props.value}
            onChange={props.onChange}>
            {props.options.map((item, index) => (
              <option key={item.value} value={item.value}> {item.text} </option>
            ))}
          </select>
        </div>
      </div>
    </React.Fragment>
  )
}
export default SelectComp;
