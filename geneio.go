// Package geneio provides interfaces for gene I/O functions.
package geneio

import (
	"fmt"
	"io"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/feat/gene"
)

const (
	maxInt = int(^uint(0) >> 1) // The maximum int value.
)

// FeaturesError describes an error that occurs when a set of features
// is associated with a gene.
type FeaturesError struct {
	Msg   string
	Gene  gene.Interface
	Feats []Feature
}

// Error returns the error message.
func (e *FeaturesError) Error() string {
	return fmt.Sprintf("%s for gene %s", e.Msg, e.Gene.Name())
}

// Reader is the common reader interface for gene.Interface.
type Reader interface {
	// Read reads a gene.Interface, returning any error that occurs during the
	// read.
	Read() (gene.Interface, error)
}

// Writer is the common writer interface for gene.Interface.
type Writer interface {
	// Write writes a gene.Interface, returning the number of bytes written
	// and any error that occurs during the write.
	Write(gene.Interface) (n int, err error)
}

// Scanner wraps a Reader to provide a convenient loop interface for reading
// genes. Successive calls to the Scan method will step through the genes of
// the provided Reader. Scanning stops unrecoverably at EOF or the first
// error.
type Scanner struct {
	r   Reader
	g   gene.Interface
	err error
}

// NewScanner returns a new Scanner to read from r.
func NewScanner(r Reader) *Scanner {
	return &Scanner{r: r}
}

// NewScannerFromFunc returns a new Scanner to read from calls to f.
func NewScannerFromFunc(f func() (gene.Interface, error)) *Scanner {
	return NewScanner(funcReader(f))
}

// Next advances the Scanner past the next gene, which will then be
// available through the Gene method. It returns false when the scan stops,
// either by reaching the end of the input or an error. After Next returns
// false, the Error method will return any error that occurred during
// scanning, except that if it was io.EOF, Error will return nil.
func (s *Scanner) Next() bool {
	if s.err != nil {
		return false
	}
	s.g, s.err = s.r.Read()
	return s.err == nil
}

// Error returns the first non-EOF error that was encountered by the Scanner.
func (s *Scanner) Error() error {
	if s.err == io.EOF {
		return nil
	}
	return s.err
}

// Gene returns the most recent gene read by a call to Next.
func (s *Scanner) Gene() gene.Interface {
	return s.g
}

// funcReader is a wrapper type that implements Reader.
type funcReader func() (gene.Interface, error)

// Read calls f and returns any gene.Interface and error returned by f.
func (f funcReader) Read() (gene.Interface, error) {
	return f()
}

// Feature is the common interface for features that can be used to construct
// a Gene.
type Feature interface {
	feat.Feature
	feat.Orienter
	// GID returns the gene ID of the Feature.
	GID() string
	// TID returns the transcript ID of the Feature.
	TID() string
	// Type returns the type of the Feature.
	// e.g. exon, start_codon, stop_codon
	Type() string
}

// FeatureReader is the common reader interface for Feature.
type FeatureReader interface {
	Read() (Feature, error)
}

// A GeneReader reads genes from a FeatureReader.
//
// It groups features with the same TID and GID into transcripts and genes
// respectively. It only considers Features of type "exon", "codon_start" and
// "codon_stop".
type GeneReader struct {
	r   FeatureReader
	blk *geneBlock
}

// NewGeneReader returns a new GeneReader that reads from r.
func NewGeneReader(r FeatureReader) *GeneReader {
	return &GeneReader{r: r}
}

// Read reads one gene from r. At EOF, it returns nil and io.EOF. When Read
// returns, r is past the last Feature incorporated in the returned gene.
func (r *GeneReader) Read() (gene.Interface, error) {
	for {
		f, err := r.r.Read()
		if err != nil {
			if err == io.EOF && r.blk != nil {
				g, err := r.blk.ToGene()
				r.blk = nil
				return g, err
			}
			return nil, err
		}
		if r.blk == nil {
			r.blk = &geneBlock{ID: f.GID(), loc: f.Location(), ori: f.Orientation()}
			r.blk.feats = append(r.blk.feats, f)
			continue
		}
		if r.blk.ID != f.GID() {
			g, err := r.blk.ToGene()
			r.blk = &geneBlock{ID: f.GID(), loc: f.Location(), ori: f.Orientation()}
			r.blk.feats = append(r.blk.feats, f)
			return g, err
		}
		r.blk.feats = append(r.blk.feats, f)
	}
}

// ReadAll reads all the remaining genes from r. A successful call returns err
// == nil, not err == io.EOF. Because ReadAll is defined to read until EOF, it
// does not treat end of file as an error to be reported. It returns a nil
// slice and an error if it encounters one.
func (r *GeneReader) ReadAll() ([]gene.Interface, error) {
	var genes []gene.Interface
	for {
		g, err := r.Read()
		if err == io.EOF {
			return genes, nil
		}
		if err != nil {
			return nil, err
		}
		genes = append(genes, g)
	}
}

// geneBlock wraps a set of Features that correspond to a single gene.
type geneBlock struct {
	ID    string
	loc   feat.Feature
	ori   feat.Orientation
	feats []Feature
}

// ToGene creates and returns a gene. It also creates the transcripts
// associated with the gene. It returns nil and an error if it encounters one.
func (geneBlk *geneBlock) ToGene() (*gene.Gene, error) {
	g := &gene.Gene{
		ID:     geneBlk.ID,
		Orient: geneBlk.ori,
		Chrom:  geneBlk.loc,
		Offset: maxInt, // is rewritten later to get the minimum start.
	}

	// Find proper gene offset and check that all features are on the same
	// location and orientation.
	for _, f := range geneBlk.feats {
		if f.Start() < g.Offset {
			g.Offset = f.Start()
		}
		if f.Orientation() != g.Orient {
			return nil, &FeaturesError{
				Gene:  g,
				Msg:   "geneio: features with varying orientation",
				Feats: geneBlk.feats,
			}
		}
		if f.Location() != g.Chrom {
			return nil, &FeaturesError{
				Gene:  g,
				Msg:   "geneio: features on varying location",
				Feats: geneBlk.feats,
			}
		}
	}

	// Build the transcripts.
	var features []feat.Feature
	var i, j int
	for i < len(geneBlk.feats) {
		tid := geneBlk.feats[i].TID()
		for j = i + 1; j < len(geneBlk.feats); j++ {
			if geneBlk.feats[j].TID() != tid {
				break
			}
		}
		tr, err := newTrancript(g, tid, geneBlk.feats[i:j])
		if err != nil {
			return nil, err
		}
		features = append(features, tr)
		i = j
	}

	// Attach the transcripts to the gene.
	if err := g.SetFeatures(features...); err != nil {
		return nil, err
	}

	return g, nil
}

// newTrancript creates and returns a new gene.Transcript from s. Returns nil
// and the error if it encounters one.
func newTrancript(g *gene.Gene, tid string, s []Feature) (gene.Transcript, error) {
	for _, f := range s {
		if f.Type() == "start_codon" {
			return newCodingTranscript(g, tid, s)
		}
	}
	return newNonCodingTranscript(g, tid, s)
}

// newCodingTranscript creates and returns a new gene.CodingTranscript from s.
// Returns nil and an error if it encounters one.
func newCodingTranscript(
	g *gene.Gene, tid string, s []Feature) (*gene.CodingTranscript, error) {

	t := &gene.CodingTranscript{ID: tid, Loc: g, Orient: feat.Forward}

	// Find offset.
	minStart := maxInt
	for _, f := range s {
		if f.Start() < minStart {
			minStart = f.Start()
		}
	}
	t.Offset = minStart - g.Offset

	// Parse Features.
	var exons []gene.Exon
	for _, f := range s {
		switch f.Type() {
		case "exon":
			e := gene.Exon{
				Transcript: t,
				Offset:     f.Start() - t.Offset - g.Offset,
				Length:     f.End() - f.Start(),
			}
			exons = append(exons, e)
		case "start_codon":
			if f.Orientation() == feat.Forward {
				t.CDSstart = f.Start() - t.Offset - g.Offset
			} else if f.Orientation() == feat.Reverse {
				t.CDSend = f.End() - t.Offset - g.Offset
			}
		case "stop_codon":
			if f.Orientation() == feat.Forward {
				t.CDSend = f.End() - t.Offset - g.Offset
			} else if f.Orientation() == feat.Reverse {
				t.CDSstart = f.Start() - t.Offset - g.Offset
			}
		}
	}

	// Set the transcript exons.
	exons = mergeExons(exons)
	if err := t.SetExons(exons...); err != nil {
		return nil, err
	}

	return t, nil
}

// newNonCodingTranscript creates and returns a new gene.NonCodingTranscript
// from s. Returns nil and an error if it encounters one.
func newNonCodingTranscript(
	g *gene.Gene, tid string, s []Feature) (*gene.NonCodingTranscript, error) {

	t := &gene.NonCodingTranscript{ID: tid, Loc: g, Orient: feat.Forward}

	// Find offset.
	minStart := maxInt
	for _, f := range s {
		if f.Start() < minStart {
			minStart = f.Start()
		}
	}
	t.Offset = minStart - g.Offset

	// Parse Features.
	var exons []gene.Exon
	for _, f := range s {
		if f.Type() == "exon" {
			e := gene.Exon{
				Transcript: t,
				Offset:     f.Start() - t.Offset - g.Offset,
				Length:     f.End() - f.Start(),
			}
			exons = append(exons, e)
		}
	}

	// Set the transcript exons.
	exons = mergeExons(exons)
	if err := t.SetExons(exons...); err != nil {
		return nil, err
	}

	return t, nil
}

// mergeExons returns a slice with touching exons concatenated into one.
// Touching exons are those that one's start equals the other's end.
func mergeExons(exons []gene.Exon) []gene.Exon {
	if len(exons) < 2 {
		return exons
	}
	var merged []gene.Exon
	e, exons := exons[0], exons[1:]
	merged = append(merged, e)
	for len(exons) > 0 {
		e, exons = exons[0], exons[1:]
		if merged[len(merged)-1].End() == e.Start() {
			merged[len(merged)-1].Length = merged[len(merged)-1].Len() + e.Len()
		} else {
			merged = append(merged, e)
		}
	}

	return merged
}
