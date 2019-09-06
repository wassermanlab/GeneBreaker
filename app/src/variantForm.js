import React from 'react';
import GeneralInfo from './generalInfo'
import VariantInfo from './variantInfo'

class VFrom extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      page: 1,
      // general state
      gene_uid: "",
      genome: "hg38",
      gene_name: "",
      chrom: "",
      sex: "XX",
      // var 1 state
      var1_type: "",
      var1_region: "",
      var1_customStart: "",
      var1_customEnd: "",
      var1_zygosity: "",
      var1_clinvar_id: "",
      var1_start: "",
      var1_end: "",
      var1_copy_change: "",
      var1_length: "",
      var1_element: "",
      var1_snv_type: "",
      var1_str_id: "",
      //var2 state
      var2: true,
      var2_type: "",
      var2_region: "",
      var2_customStart: "",
      var2_customEnd: "",
      var2_zygosity: "",
      var2_clinvar_id: "",
      var2_start: "",
      var2_end: "",
      var2_copy_change: "",
      var2_length: "",
      var2_element: "",
      var2_snv_type: "",
      var2_str_id: "",
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
      const index = target.selectedIndex;
      const option = target.childNodes[index];
      const chrom = option.getAttribute('chrom');
      this.setState({
        gene_uid: value,
        chrom: chrom
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
        {/* //page, chrom, sex, gene_uid, handleInputChange, next, back
        //var, var2, type, region, zygosity, clinvar_id, start, end, copy_change, length, element, snv_type, motif */}
        <VariantInfo 
          page={this.state.page} 
          chrom={this.state.chrom} sex={this.state.sex} gene_uid={this.state.gene_uid} 
          handleInputChange={this.handleInputChange} next={this.next}  back={this.back}
          var={1} type={this.state.var1_type} region={this.state.var1_region} zygosity={this.state.var1_zygosity}
          customStart={this.state.var1_customStart} customEnd={this.state.var1_customEnd}
          start={this.state.var1_start} end={this.state.var1_end} length={this.state.var1_length}
          copy_change={this.state.var1_copy_change} clinvar_id={this.state.var1_clinvar_id} 
          element={this.state.var1_element} snv_type={this.state.var1_snv_type} str_id={this.state.var1_str_id}
        />
        {/* var2 */}
         {/* <VariantInfo 
          page={this.state.page} 
          chrom={this.state.chrom} sex={this.state.sex} gene_uid={this.state.gene_uid} 
          handleInputChange={this.handleInputChange} next={this.next}  back={this.back}
          var={2} var2={this.state.var2}
          type={this.state.var2_type} region={this.state.var2_region} zygosity={this.state.var2_zygosity}
          customStart={this.state.var2_customStart} customStart={this.state.var2_customStart}
          start={this.state.var2_start} end={this.state.var2_end} length={this.state.var2_length}
          copy_change={this.state.var2_copy_change} clinvar_id={this.state.var2_clinvar_id} 
          element={this.state.var2_element} snv_type={this.state.var2_snv_type} motif={this.state.var2_motif}
        /> */}
      </form>
    );
  }
}

export default VFrom;
