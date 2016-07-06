// Package gff reads genes from a GFF v2 file.
//
// It requires that GFF entries, particularly entries of type exon,
// start_codon and stop_codon are grouped into genes and transcripts based on
// GFF tags. Entries are expected to be sorted by gene and transcript grouping
// tags as in the example below.
//
//  Y	.	exon	10	20	0	-	.	gene_id A; transcript_id A1
//  Y	.	exon	50	90	0	-	.	gene_id A; transcript_id A1;
//  Y	.	stop_codon	60	62	0	-	.	gene_id A; transcript_id A1
//  Y	.	start_codon	71	73	0	-	.	gene_id A; transcript_id A1
//  Y	.	exon	10	100	0	-	.	gene_id A; transcript_id A2
//  Y	.	stop_codon	60	62	0	-	.	gene_id A; transcript_id A2
//  Y	.	start_codon	71	73	0	-	.	gene_id A; transcript_id A2
//
package gff

import (
	"errors"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/feat/gene"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/go-bio/geneio"
)

// A Reader reads genes from a GFF v2 file.
//
// It groups features with the same transcript and gene group tag values into
// transcripts and genes respectively. Group tags can be changed by SetGeneTag
// and SetTranscriptTag to customize the details before the first call to Read
// or ReadAll.
type Reader struct {
	r         *geneio.GeneReader
	fr        *featureReader
	afterRead bool
}

// NewReader returns a new Reader that reads from r. It sets transcript and
// gene group tag to "transcript_id" and "gene_id" respectively.
func NewReader(r featio.Reader) *Reader {
	fr := &featureReader{
		r:             r,
		GeneTag:       "gene_id",
		TranscriptTag: "transcript_id",
	}
	return &Reader{r: geneio.NewGeneReader(fr), fr: fr}
}

// Read reads one gene from r. At EOF, it returns nil and io.EOF. When Read
// returns, r is past the last feature incorporated in the returned gene.
func (r *Reader) Read() (gene.Interface, error) {
	r.afterRead = true
	return r.r.Read()
}

// ReadAll reads all the remaining genes from r. A successful call returns err
// == nil, not err == io.EOF. Because ReadAll is defined to read until EOF, it
// does not treat end of file as an error to be reported. It returns a nil
// slice and an error if it encounters one.
func (r *Reader) ReadAll() ([]gene.Interface, error) {
	r.afterRead = true
	return r.r.ReadAll()
}

// SetGeneTag sets the gene group tag. It is set to ('gene_id') by NewReader.
// Tag can only be changed before the first call to Read or ReadAll.
func (r *Reader) SetGeneTag(tag string) error {
	if r.afterRead {
		return errors.New("gff: cannot set GeneTag after first call to Read")
	}
	r.fr.GeneTag = tag
	return nil
}

// SetTranscriptTag sets the transcript group tag. It is set to
// ('transcript_id') by NewReader. Tag can only be changed before the first
// call to Read or ReadAll.
func (r *Reader) SetTranscriptTag(tag string) error {
	if r.afterRead {
		return errors.New("gff: cannot set GeneTag after first call to Read")
	}
	r.fr.TranscriptTag = tag
	return nil
}

// featureReader is an implementation of geneio.FeatureReader.
type featureReader struct {
	r                      featio.Reader
	GeneTag, TranscriptTag string
}

// Read implements geneio.FeatureReader.
func (r *featureReader) Read() (geneio.Feature, error) {
	f, err := r.r.Read()
	if err != nil {
		return nil, err
	}
	return r.NewFeature(f)
}

// NewFeature converts f to *feature and returns it.
func (r *featureReader) NewFeature(f feat.Feature) (*feature, error) {
	gf, ok := f.(*gff.Feature)
	if !ok {
		return nil, errors.New("gff: feature is not a GFF feature")
	}
	gid := gf.FeatAttributes.Get(r.GeneTag)
	if gid == "" {
		return nil, errors.New("gff: empty grouping " + r.GeneTag + " field")
	}
	tid := gf.FeatAttributes.Get(r.TranscriptTag)
	if tid == "" {
		return nil, errors.New("gff: empty grouping " + r.TranscriptTag + " field")
	}

	return &feature{
		Feature: f,
		fgid:    gid,
		ftid:    tid,
		ftype:   gf.Feature,
		ori:     feat.Orientation(gf.FeatStrand),
	}, nil
}

// feature is an implementation of geneio.Feature.
type feature struct {
	feat.Feature
	ori               feat.Orientation
	fgid, ftid, ftype string
}

func (f *feature) GID() string                   { return f.fgid }
func (f *feature) TID() string                   { return f.ftid }
func (f *feature) Type() string                  { return f.ftype }
func (f *feature) Orientation() feat.Orientation { return f.ori }
