import React from 'react';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import Form from 'react-bootstrap/Form';
import InputGroup from 'react-bootstrap/InputGroup';
import Button from 'react-bootstrap/Button';

class General extends React.Component {
    render() {
        if (this.props.currentStep !== 0) { // Prop: The current step
            return null
        }
        // The markup for the Step 1/2 UI
        return (
            <div>
                <h2>General Information</h2>
                <Form.Group as={Row} controlId="genome">
                    <Form.Label column sm={2}>Genome</Form.Label>
                    <Col sm={10}>
                        <Form.Control as="select">
                            <option>hg38</option>
                            <option>hg19</option>
                        </Form.Control>
                    </Col>
                </Form.Group>
                <Form.Group as={Row} controlId="sex">
                    <Form.Label column sm={2}>Proband Sex</Form.Label>
                    <Col sm={10}>
                        <Form.Control as="select">
                            <option>XX</option>
                            <option>XY</option>
                        </Form.Control>
                    </Col>
                </Form.Group>
                <Form.Group as={Row} controlId="geneName">
                    <Form.Label column sm={2}>
                        Gene Name
                    </Form.Label>
                    <Col sm={10}>
                        <InputGroup>
                            <Form.Control type="text" />
                            <InputGroup.Append>
                                <Button variant="outline-secondary">Get Transcripts</Button>
                            </InputGroup.Append>
                        </InputGroup>
                    </Col>
                </Form.Group>
                <Form.Group as={Row} controlId="transcripts">
                    <Form.Label column sm={2}>Transcripts</Form.Label>
                    <Col sm={10}>
                        <Form.Control as="select">
                        </Form.Control>
                    </Col>
                </Form.Group>
                <Row>
                <Col sm={11}></Col>
                <Col sm={1}>
                    <Button  variant="primary" >Next</Button></Col>
                </Row>
            </div>

        )
    }
}

export default General;


