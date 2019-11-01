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
          </tr>
        </thead>
        <tbody>
          <tr>
            <th scope="row">1</th>
            <td>XX</td>
            <td><input type="checkbox" /></td>
            <td><input type="checkbox" /></td>
            <td><input type="checkbox" /></td>
          </tr>
        </tbody>
      </table>
    </React.Fragment >
  )
}


export default FamilyInfo;