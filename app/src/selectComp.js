import React from 'react';

function SelectComp(props) {
  return (
    <React.Fragment>
      <div className="form-group">
        <label>{props.title}</label>
        <select className="form-control"
          name={props.name}
          value={props.value}
          onChange={props.onChange}>
          {props.options.map((item, index) => (
              <option key={item.value} value={item.value}> {item.text} </option>
          ))}
      </select>
      </div>
    </React.Fragment>
  )
}
export default SelectComp;


// <SelectComp
//           title={"Genome"}
//           name={"genome"}
//           value={this.props.genome}
//           onChange={this.props.handleInputChange}
//           options={[{value: "hg38", text: "hg38"},{value: "hg19", text: "hg19"}]}
//           />