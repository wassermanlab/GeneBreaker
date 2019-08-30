import React from 'react';
import Form from 'react-bootstrap/Form';
import FormControl from 'react-bootstrap/FormControl';
import Button from 'react-bootstrap/Button';
import InputGroup from 'react-bootstrap/InputGroup';

class VFrom2 extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      gene_uid: "",
      genome: "hg38",
      gene_name: "",
      chr: "",
      sex: "XX",
      var1: {
        type: "",
        region: "",
        zygosity: "",
        impact: {
          clinvar: {
            id: ""
          },
          cnv: {
            start: "",
            end: "",
            copy_change: ""
          },
          indel: {
            length: "",
            start: ""
          },
          mei: {
            element: "",
            start: ""
          },
          snv: {
            snv_type: "",
            start: ""
          },
          str: {
            start: "",
            end: "",
            motif: "",
            length: ""
          }
        }
      },
      var2: {
        type: "",
        region: "",
        zygosity: "",
        impact: {
          clinvar: {
            id: ""
          },
          cnv: {
            start: "",
            end: "",
            copy_change: ""
          },
          indel: {
            length: "",
            start: ""
          },
          mei: {
            element: "",
            start: ""
          },
          snv: {
            snv_type: "",
            start: ""
          },
          str: {
            start: "",
            end: "",
            motif: "",
            length: ""
          }
        }
      }
    };

    this.handleInputChange = this.handleInputChange.bind(this);
  }

  handleInputChange(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;

    this.setState({
      [name]: value
    },
      () => console.log(this.state));
  }

  render() {
    return (
      <Form>
        {/* general: genome, sex, gene_uid, chr, transcript */}
        <React.Fragment>
          <Form.Group controlId="exampleForm.ControlSelect1">
            <Form.Label>Genome</Form.Label>
            <Form.Control
              as="select"
              type="select"
              name="genome"
              value={this.state.genome}
              onChange={this.handleInputChange} >
              <option value="hg38">hg38</option>
              <option value="hg19">hg19</option>
            </Form.Control>
          </Form.Group>

          <Form.Group controlId="exampleForm.ControlSelect1">
            <Form.Label>Sex</Form.Label>
            <Form.Control
              as="select"
              type="select"
              name="sex"
              value={this.state.sex}
              onChange={this.handleInputChange} >
              <option value="XX">XX</option>
              <option value="XY">XY</option>
            </Form.Control>
          </Form.Group>

          <label>Gene Name</label>
          <InputGroup className="mb-3">
            <FormControl name="gene_name" value={this.state.gene_name} onChange={this.handleInputChange} />
            <InputGroup.Append>
              <Button variant="outline-secondary">Get Transcripts</Button>
            </InputGroup.Append>
          </InputGroup>
          {/* [
  {
    "chrom": "chr17", 
    "end": 72126413, 
    "id": "genes_69540", 
    "qualifiers": {
      "cdsEnd": 72124387, 
      "cdsStart": 72121391, 
      "exonEnds": "72121822,72122972,72126413,", 
      "exonStarts": "72121019,72122718,72123542,", 
      "name": "NM_000346", 
      "name2": "SOX9", 
      "source": "refGene", 
      "uid": 69540
    }, 
    "score": 0, 
    "start": 72121019, 
    "strand": 1
  }
] */}


        </React.Fragment>
        {/* var1 */}
        {/* var2 */}
      </Form>
    );
  }
}

export default VFrom2;
