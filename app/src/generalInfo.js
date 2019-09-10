import React from 'react';
import NavButtons from './navButtons'
import SelectComp from './selectComp'

// todo: move to separate component 
function GeneralErrors(props) {
  const select_error = "No transcript is selected.";
  const sex_error = "Transcript on Y chromosome cannot be selected with XX proband."
  if (props.transcript_select) {
    return <div className="alert alert-danger" role="alert">{select_error}</div>
  }
  else if (props.transcript_sex) {
    return <div className="alert alert-danger" role="alert">{sex_error}</div>
  }
  return null;
}

class GeneralInfo extends React.Component {
  //page, chrom, genome, sex, gene_name, gene_uid, handleInputChange, next, back

  constructor(props) {
    super(props);
    this.state = {
      transcript_list: [],
      no_transcript: false,
      transcript_select: false,
      transcript_sex: false
    };
    this.getTranscripts = this.getTranscripts.bind(this);
    this.customNext = this.customNext.bind(this);

  }

  async getTranscripts(event) {
    const response = await fetch('http://127.0.0.1:5001/get_transcripts/' + this.props.genome + '/' + this.props.gene_name);
    const json = await response.json();
    if (!json) {
      this.setState({ no_transcript: true, transcript_list: [] })
    }
    this.setState({ no_transcript: false, transcript_list: json })
  }

  // error check for the whole page!
  customNext() {
    if (this.props.gene_uid === "") {
      this.setState({ transcript_select: true })
    } else if (this.props.chrom === "chrY" && this.props.sex === "XX") {
      this.setState({ transcript_sex: true })
    } else {
      this.props.next()
    }
  }

  render() {
    if (this.props.page !== 1) {
      return null;
    }
    return (

      // {/* general: genome, sex, gene_uid, chr, transcript */}
      <React.Fragment>
        {/* genome */}
        <SelectComp
          title={"Genome"}
          name={"genome"}
          value={this.props.genome}
          onChange={this.props.handleInputChange}
          options={[{value: "hg38", text: "hg38"},{value: "hg19", text: "hg19"}]}
          />
        {/* sex */}
        <SelectComp
          title={"Sex"}
          name={"sez"}
          value={this.props.sex}
          onChange={this.props.handleInputChange}
          options={[{value: "XX", text: "XX"},{value: "XY", text: "XY"}]}
          />
        {/* gene name/symbol */}
        <label>Gene symbol</label>
        <div className="input-group mb-3">
          <input type="text" className="form-control" name="gene_name" value={this.props.gene_name} onChange={this.props.handleInputChange} />
          <div className="input-group-append">
            <button className="btn btn-outline-secondary" type="button" onClick={this.getTranscripts}>Get Transcripts</button>
          </div>
        </div>
        {/* gene/name errors */}
        {this.state.no_transcript ?
          <div className="alert alert-danger" role="alert">No transcripts found by the requested name.</div> : null
        }
        {/* transcript list */}
        <div className="form-group">
          <label>Transcript</label>
          <select className="form-control"
            name="gene_uid"
            value={this.props.gene_uid}
            onChange={this.props.handleInputChange}
            size="5">
            <option key="0" value=""></option>
            {this.state.transcript_list.length > 0 ? (this.state.transcript_list.map((item, index) => (
              <option key={item.qualifiers.uid} value={item.qualifiers.uid} chrom={item.chrom}> {item.qualifiers.name} </option>
            ))) : null}
          </select>
        </div>
        <GeneralErrors
          transcript_select={this.state.transcript_select}
          transcript_sex={this.state.transcript_sex}
        />
        <NavButtons
          next={this.customNext}
          back={this.props.back}
          page={this.props.page} />
      </React.Fragment>
    );
  }


}

export default GeneralInfo;
