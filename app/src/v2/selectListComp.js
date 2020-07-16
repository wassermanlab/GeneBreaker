import React, { useState } from 'react';

function SelectList(props) {
  const [populatingText, setPopulatingText] = useState('')
  const [list, setList] = useState([]);

  async function populateBox() {
    // set the text 
    setPopulatingText('Populating with ' + props.populatingText)
    try {
      const response = await fetch(props.url)
      const json = await response.json();
      setList(json);
      if (json.length > 0) {
        setPopulatingText('FILLED')
      }
      else {
        setPopulatingText('No results found with with ' + props.populatingText)
      }
    } catch (error) {
      console.error('Error:', error)
      setPopulatingText('Error upon populating.')
    }
  }

  function option(item, index) {
    if (props.type === "transcript") {
      return <option key={index}
        value={item.qualifiers.uid}
        start={item.start + 1}
        end={item.end}
        chrom={item.chrom}>
        {item.qualifiers.accession_number} {item.qualifiers.source === "ncbiRefSeqSelect"? "*" :""} </option>
    } else if (props.type === "clingen") {
      return (<option key={index} value={item.qualifiers.uid} start={item.start + 1} end={item.end} copy_change={item.qualifiers.copy_number_change}>
        Start: {item.start + 1} &emsp;&emsp; End: {item.end} &emsp;&emsp; Clinical assertion: {item.qualifiers.clinical_assertion} &emsp;&emsp; Copy change: {item.qualifiers.copy_number_change} </option>
      )
    } else if (props.type === 'str') {
      return (<option key={index} value={item.qualifiers.uid}>
        Start: {item.start + 1} &emsp;&emsp; End: {item.end} &emsp;&emsp; Motif: {item.qualifiers.motif} &emsp;&emsp; Pathogenicity: {item.qualifiers.pathogenicity} </option>)
    } else if (props.type === 'clinvar') {
      return (<option key={index} value={item.qualifiers.clinvarID}>
        ClinVarID: {item.qualifiers.clinvarID}  &emsp;&emsp; Start: {item.start + 1} &emsp;&emsp; Ref: {item.qualifiers.ref} &emsp;&emsp; Alt: {item.qualifiers.alt} &emsp;&emsp; Clinical Significance: {item.qualifiers.CLNSIG} </option>
      )
    }
  }

  return (
    <React.Fragment>
      <div className="form-group row">
        <label className="col-sm-2 col-form-label" >{props.title}</label>
        <div className="col-sm-10">

          {/* populate button */}
          <button onClick={populateBox} type="button"
            className="btn btn-light btn-sm"
            style={{ marginBottom: "10px" }}>
            Fetch {props.title}
          </button>

          <br />

          {/* populating message */}
          {populatingText !== "FILLED" &&
            <small className="form-text text-muted">
              {populatingText}
            </small>}

          {/* list */}

          {populatingText === "FILLED" &&
            <select className="form-control"
              name={props.name}
              value={props.value}
              onChange={props.onChange}
              size="5">
              <option key="0" value=""></option>
              {/* map of options*/}
              {list.map((item, index) => (
                option(item, index)
              ))}
            </select>
          }
          {props.type === "transcript" && populatingText === "FILLED" && 
          <small className="form-text text-muted">
          Asterisk indicates the representative transcript as chosen by the RefSeqSelect 
          dataset, for more info please refer to the documentation on 
          <a href="https://www.ncbi.nlm.nih.gov/refseq/refseq_select/"> NCBI RefSeq Select.</a> 
        </small>
          }

        </div>
      </div>
    </React.Fragment>
  )

}

export default SelectList;