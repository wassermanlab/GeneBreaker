import React from 'react';
import GeneralInfo from './generalInfo'
import NavButtons from './navButtons'

class VFrom extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      page: 1,
      clinvar_list: [],
      clingen_list: [],
      str_list: [],
      gene_uid: "",
      genome: "hg38",
      gene_name: "",
      chrom: "",
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
    this.next = this.next.bind(this)
    this.back = this.back.bind(this)
  }

  // sets page to page+1
  next(){
    let currentPage = this.state.page;

    currentPage = currentPage + 1; 
    this.setState({page: currentPage});
  }

  // sets page to page-1
  back(){
    let currentPage = this.state.page;

    currentPage = currentPage - 1; 
    this.setState({page: currentPage});
  }

  handleInputChange(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;

    //custom state change to deal with chrom
    if (name === "gene_uid" && value !== "") {
      const parts = value.split("_");
      this.setState({
        gene_uid: parts[1],
        chrom: parts[0]
      },
        () => console.log(this.state));
      return null;
    }
    // state change for everything else 
    this.setState({
      [name]: value
    },
      () => console.log(this.state));
  }
  render() {
    return (
      <form>
        {/* general: genome, sex, gene_uid, chr, transcript */}
        <GeneralInfo 
          page={this.state.page}
          chrom={this.state.chrom}
          genome={this.state.genome}
          sex={this.state.sex}
          gene_name={this.state.gene_name}
          gene_uid={this.state.gene_uid}
          handleInputChange={this.handleInputChange}
          next={this.next} 
          back={this.back}
        />
        {/* var1 */}
        {/* var2 */}
      </form>
    );
  }
}

export default VFrom;
