import React from 'react';
import Container from 'react-bootstrap/Container';
import Form from 'react-bootstrap/Form';
import General from './variantsGeneral'
import Family from './variantsFamily'
import Variant from './variantsVariant'

// class DesignVariants extends React.Component {

//   constructor(props) {
//     super(props);
//     this.state = {
//       currentStep: 0,         // panel options: 0=generalInfo, 1=var1, 2=var2, 3=family
//       genome: "hg38",
//       sex: "XX",
//       chromosome: "",
//       gene_uid: "",
//       var1: {
//         type: "",
//         region: "",
//         zygosity: "",
//         clinvar_id: "",
//         clinvar_start: "",
//         clinvar_ref: "",
//         clinvar_alt: "",
//         cnv_chrom: "",
//         cnv_start: "",
//         cnv_end: "",
//         cnv_copy_change: "",
//         indel_amount: "",
//         indel_start: "",
//         mei_element: "",
//         mei_start: "",
//         str_chrom: "",
//         str_start: "",
//         str_end: "",
//         str_length: "",
//         str_motif: "",
//         snv_type: "",
//         snv_start: ""
//       },
//       var2: {
//         type: "",
//         region: "",
//         zygosity: "",
//         clinvar_id: "",
//         clinvar_start: "",
//         clinvar_ref: "",
//         clinvar_alt: "",
//         cnv_chrom: "",
//         cnv_start: "",
//         cnv_end: "",
//         cnv_copy_change: "",
//         indel_amount: "",
//         indel_start: "",
//         mei_element: "",
//         mei_start: "",
//         str_chrom: "",
//         str_start: "",
//         str_end: "",
//         str_length: "",
//         str_motif: "",
//         snv_type: "",
//         snv_start: ""
//       }
//     };
//     // Bind the submission to handleChange() 
//     this.handleChange = this.handleChange.bind(this)
//   }

//   render() {
//     return (

//       // name="numberOfGuests"
//       //         type="number"
//       //         value={this.state.numberOfGuests}
//       //         onChange={this.handleInputChange}
//       <React.Fragment>
//         <Container>
//           <Form onSubmit={this.handleSubmit}>
//             {/* general info */}
//             <General currentStep={this.state.currentStep} />
//             {/* variant 1 */}
//             <Variant currentStep={this.state.currentStep} />
//             {/* variant 2 */}
//             <Variant currentStep={this.state.currentStep} />
//             {/* family info */}
//             <Family currentStep={this.state.currentStep} />

//           </Form>
//         </Container>
//       </React.Fragment>
//     );
//   }
// }

// export default DesignVariants;


class DesignVariants extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      currentStep: 0,
      genome: 'hg38',
      sex: 'XX',
      gene_name: "",
      gene_uid: "",
      transcripts: [
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
      ]
    };

    this.handleInputChange = this.handleInputChange.bind(this);
    this.getTranscripts = this.getTranscripts.bind(this);
  }

  getTranscripts(event) {
    const axios = require('axios');
    axios.get('http://127.0.0.1:5001/get_transcripts/' + this.state.genome + '/' + this.state.gene_name)
      .then(function (response) {
        // populate transcripts
        console.log(response.data);;
      })
      .catch(function (error) {
        // print error
        return (error);
      })
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
      <React.Fragment>
        <Container>
          <Form onSubmit={this.handleSubmit}>
            {/* general info */}
            <General
              currentStep={this.state.currentStep}
              handleInputChange={this.handleInputChange}
              getTranscripts={this.getTranscripts}
              genome={this.state.genome}
              sex={this.state.sex}
              gene_name={this.state.gene_name}
              gene_uid={this.state.gene_uid}
              transcripts={this.state.transcripts}
            />
            {/* variant 1 */}
            <Variant currentStep={this.state.currentStep} />
            {/* variant 2 */}
            <Variant currentStep={this.state.currentStep} />
            {/* family info */}
            <Family currentStep={this.state.currentStep} />
          </Form>
        </Container>
      </React.Fragment>
    );
  }
}

export default DesignVariants;