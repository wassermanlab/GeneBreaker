import React from 'react';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import Form from 'react-bootstrap/Form';
import InputGroup from 'react-bootstrap/InputGroup';
import Button from 'react-bootstrap/Button';
class Variant extends React.Component {
  render() {
    if (this.props.currentStep > 2 || this.props.currentStep < 1) { // Prop: The current step
      return null
    }
    // The markup for the Step 1/2 UI
    return (
      <div>
        <h2>Variant {this.props.currentStep}</h2>
        <Form.Group as={Row} controlId="type">
          <Form.Label column sm={2}>Type</Form.Label>
          <Col sm={10}>
            <Form.Control as="select">
              <option>Select</option>
              <option>ClinVar Variant</option>
              <option>ClinGen Copy Number Variant</option>
              <option>Copy Number Variant</option>
              <option>Indel</option>
              <option>Mobile Element Insertion</option>
              <option>Single Nucleotide Variant</option>
              <option>Short Tandem Repeat</option>
            </Form.Control>
          </Col>
        </Form.Group>
        <Form.Group as={Row} controlId="region">
          <Form.Label column sm={2}>Region</Form.Label>
          <Col sm={10}>
            <Form.Control as="select">
              <option>Select</option>
              <option>Genic</option>
              <option>Coding</option>
              <option>Intronic</option>
              <option>Untranslated Region</option>
              <option>Custom</option>
            </Form.Control>
          </Col>
        </Form.Group>
        <Form.Group as={Row} controlId="zygosity">
          <Form.Label column sm={2}>Zygosity</Form.Label>
          <Col sm={10}>
            <Form.Control as="select">
              <option>Select</option>
              <option>Homozygous</option>
              <option>Heterozygous</option>
              <option>Hemizygous</option>
            </Form.Control>
          </Col>
        </Form.Group>
        <Row>
          <Col sm={10}></Col>
          <Col sm={1}>
            <Button variant="secondary">Back</Button>
          </Col>
          <Col sm={1}>
            <Button variant="primary">Next</Button>
          </Col>
        </Row>
      </div>
    )
  }
}


export default Variant;


