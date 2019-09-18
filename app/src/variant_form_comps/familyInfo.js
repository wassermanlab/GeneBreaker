import React from 'react';

function FInfo(props) {
    // element, start
    if (props.page !== 4) {
        return null;
    }
    return (<React.Fragment>
        <h2>Variants</h2>
        <table className="table">
            <thead>
                <tr>
                    <th scope="col">Variant</th>
                    <th scope="col">CHROM</th>
                    <th scope="col">POS</th>
                    <th scope="col">ID</th>
                    <th scope="col">REF</th>
                    <th scope="col">ALT</th>
                    <th scope="col">QUAL</th>
                    <th scope="col">FILTER</th>
                    <th scope="col">INFO</th>
                    <th scope="col">FORMAT</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <th scope="row">1</th>
                    <td>{props.vars.var1.chrom}</td>
                    <td>{props.vars.var1.pos}</td>
                    <td>{props.vars.var1.id}</td>
                    <td>{props.vars.var1.ref}</td>
                    <td>{props.vars.var1.alt}</td>
                    <td>{props.vars.var1.qual}</td>
                    <td>{props.vars.var1.filter}</td>
                    <td>{props.vars.var1.info}</td>
                    <td>{props.vars.var1.format}</td>
                </tr>
                {props.vars.var1.proband === "1/1" ?
                    <tr>
                        <th scope="row">2</th>
                        <td colSpan="9">Homozygous with variant 1</td>
                    </tr> : null}
                {(props.vars.var1.proband === "0/1" && props.vars.var2 !== "") ?
                    <tr>
                        <th scope="row">2</th>
                        <td colSpan="9">No variant 2</td>
                    </tr> : null}
                {props.vars.var2 !== "" ? (
                    <tr>
                        <th scope="row">1</th>
                        <td>{props.vars.var2.chrom}</td>
                        <td>{props.vars.var2.pos}</td>
                        <td>{props.vars.var2.id}</td>
                        <td>{props.vars.var2.ref}</td>
                        <td>{props.vars.var2.alt}</td>
                        <td>{props.vars.var2.qual}</td>
                        <td>{props.vars.var2.filter}</td>
                        <td>{props.vars.var2.info}</td>
                        <td>{props.vars.var2.format}</td>
                    </tr>) : null}
            </tbody>
        </table>
        <h2>Family</h2>
        {/* id -- sex -- var1 -- var2 -- affected */}
        <table className="table">
            <thead>
                <tr>
                    <th scope="col">Family Member</th>
                    <th scope="col">Sex</th>
                    <th scope="col">Variant 1</th>
                    <th scope="col">Variant 2</th>
                    <th scope="col">Affected</th>
                    <th scope="col">Remove</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <th scope="row">Proband</th>
                    <td>{props.sex}</td>
                    <td><input type="checkbox" checked disabled /></td>
                    <td>
                        {((props.vars.var1.proband === "0/1" && props.vars.var2 !== "") || 
                        (["chrX", "chrY"].includes(props.vars.var1.proband.chrom) && props.sex === "XY" )) ? 
                        <input type="checkbox" disabled /> : <input type="checkbox" checked disabled />}
                    </td>
                    <td><input type="checkbox" checked disabled /></td>
                    <td><button type="button" className="btn btn-primary" disabled>Remove</button></td>
                </tr>
                {Object.keys(props.family).map((item, index) => (
                    <tr key={item}>
                        <th scope="row">{item}</th>
                        <td>{props.family[item].sex}</td>
                        <td><input type="checkbox" onChange={props.onChange} name={"var1_" + item} /></td>
                        <td><input type="checkbox" onChange={props.onChange} name={"var2_" + item} /></td>
                        <td><input type="checkbox" onChange={props.onChange} name={"affected_" + item} /></td>
                        <td><button type="button" className="btn btn-primary" onClick={props.onRemove} value={item} >Remove</button></td>
                    </tr>
                ))}
            </tbody>
        </table>
        <div className="dropdown">
            <button className="btn btn-primary dropdown-toggle float-right" type="button" data-toggle="dropdown">
                Add family Member</button>
            <div className="dropdown-menu">
                <button className="dropdown-item" type="button" onClick={props.onAdd} value={"m"} >Mother</button>
                <button className="dropdown-item" type="button" onClick={props.onAdd} value={"f"} >Father</button>
                <button className="dropdown-item" type="button" onClick={props.onAdd} value={"s"} >Sister</button>
                <button className="dropdown-item" type="button" onClick={props.onAdd} value={"b"} >Brother</button>
            </div>
        </div>
        <div className="dropdown">
            <button className="btn btn-primary dropdown-toggle float-right" type="button" data-toggle="dropdown">
                Download output</button>
            <div className="dropdown-menu">
                <button className="dropdown-item" type="button" onClick={props.downloadFile} value={"vcf"} >VCF</button>
                <button className="dropdown-item" type="button" onClick={props.downloadFile} value={"ped"} >PED</button>
            </div>
        </div>
    </React.Fragment>)
}
export default FInfo;
