import React from 'react';
// import InputComp from './inputComp'
// import SelectComp from './selectComp'
// import SelectList from './selectListComp'
// import { host } from '../host'

function FamilyInfo(props) {
  if (props.page !== 4) {
    return null;
  }
  let row2;
  if (props.vars.var2 === "") {
    if (props.vars.var1.proband === "1/1") {
      row2 = (
        <tr>
          <th scope="row">2</th>
          <td colSpan="4">Same as variant 1</td>
        </tr>)
    } else {
      row2 = (
        <tr>
          <th scope="row">2</th>
          <td colSpan="4">None</td>
        </tr>)
    }
  }
  else {
    row2 = (
      <tr>
        <th scope="row">2</th>
        <td>{props.vars.var2.chrom}</td>
        <td>{props.vars.var2.pos}</td>
        <td>{props.vars.var2.ref}</td>
        <td>{props.vars.var2.alt}</td>
      </tr>)
  }
  return (
    <React.Fragment>
      <h1>Variants</h1>
      <table className="table">
        <thead>
          <tr>
            <th scope="col">Variant</th>
            <th scope="col">Chrom</th>
            <th scope="col">Pos</th>
            <th scope="col">Ref</th>
            <th scope="col">Alt</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th scope="row">1</th>
            <td>{props.vars.var1.chrom}</td>
            <td>{props.vars.var1.pos}</td>
            <td>{props.vars.var1.ref}</td>
            <td>{props.vars.var1.alt}</td>
          </tr>
          {row2}
        </tbody>
      </table>
      <h1>Family Info</h1>
      <table className="table">
        <thead>
          <tr>
            <th scope="col">ID</th>
            <th scope="col">Sex</th>
            <th scope="col">Var1</th>
            <th scope="col">Var2</th>
            <th scope="col">Affected</th>
            <th scope="col"></th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th scope="row">proband</th>
            <td>{props.family.proband.sex}</td>
            <td><input type="checkbox" checked disabled /></td>
            <td>
              {((props.vars.var1.proband === "0/1" && props.vars.var2 === "") ||
                (["chrX", "chrY"].includes(props.vars.var1.chrom) && props.family.proband.sex === "XY")) ?
                <input type="checkbox" disabled /> : <input type="checkbox" checked disabled />}
            </td>
            <td><input type="checkbox" checked disabled /></td>
            <td><button type="button" className="btn btn-secondary" disabled>Remove</button></td>
          </tr>
          {Object.keys(props.family).map((item, index) => (
            (item !== "proband" ? <tr key={item}>
            <th scope="row">{item}</th>
            <td>{props.family[item].sex}</td>
            <td><input type="checkbox" onChange={props.onChange} name={"var1_" + item} /></td>
            <td><input type="checkbox" onChange={props.onChange} name={"var2_" + item} /></td>
            <td><input type="checkbox" onChange={props.onChange} name={"affected_" + item} /></td>
            <td><button type="button" className="btn btn-secondary" onClick={props.onRemove} value={item} >Remove</button></td>
          </tr>: null)
            
          ))}
        </tbody>
      </table>
      <div>
        <div className="dropdown">
          <button className="btn btn-primary float-left dropdown-toggle" type="button" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
            Add Family Member
          </button>
          <div className="dropdown-menu" aria-labelledby="dropdownMenuButton">
            <button className="dropdown-item" type="button" value={"m"} >Mother</button>
            <button className="dropdown-item" type="button" value={"f"} >Father</button>
            <button className="dropdown-item" type="button" value={"s"} >Sister</button>
            <button className="dropdown-item" type="button" value={"b"} >Brother</button>
          </div>
        </div>
      </div>
    </React.Fragment >
  )
}


export default FamilyInfo;