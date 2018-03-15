package geneio

import (
	"errors"
	"io"
	"reflect"
	"sort"
	"strings"
	"testing"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/feat/gene"
	"github.com/biogo/biogo/io/featio/gff"
)

// Assert that interfaces are satisfied.
var (
	_ Reader        = (*GeneReader)(nil)
	_ FeatureReader = (*FeatureReaderImpl)(nil)
)

// Define this type to enable sorting and make the tests stable.
type Genes []gene.Interface

func (s Genes) Len() int {
	return len(s)
}

func (s Genes) Less(i, j int) bool {
	return s[i].Name() < s[j].Name()
}

func (s Genes) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

// FeatureReaderImpl is a simple implementation of FeatureReader.
type FeatureReaderImpl struct {
	r        *gff.Reader
	GID, TID string
}

func (r *FeatureReaderImpl) Read() (Feature, error) {
	f, err := r.r.Read()
	if err != nil {
		return nil, err
	}
	return r.NewFeatureImpl(f)
}

func (r *FeatureReaderImpl) NewFeatureImpl(f feat.Feature) (*FeatureImpl, error) {
	gf, ok := f.(*gff.Feature)
	if !ok {
		return nil, errors.New("gff: feature is not a GFF feature")
	}
	gid := gf.FeatAttributes.Get(r.GID)
	if gid == "" {
		return nil, errors.New("gff: empty grouping " + r.GID + " field")
	}
	tid := gf.FeatAttributes.Get(r.TID)
	if tid == "" {
		return nil, errors.New("gff: empty grouping " + r.TID + " field")
	}
	return &FeatureImpl{
		Feature: f,
		fgid:    gid,
		ftid:    tid,
		ftype:   gf.Feature,
		ori:     feat.Orientation(gf.FeatStrand),
	}, nil
}

// FeatureImpl is a simple implementation of Feature.
type FeatureImpl struct {
	feat.Feature
	ori               feat.Orientation
	fgid, ftid, ftype string
}

func (f *FeatureImpl) GID() string                   { return f.fgid }
func (f *FeatureImpl) TID() string                   { return f.ftid }
func (f *FeatureImpl) Type() string                  { return f.ftype }
func (f *FeatureImpl) Orientation() feat.Orientation { return f.ori }

// geneGenerator is the simplest gene generator. It just returns the gene that
// it gets.
func geneGenerator(g gene.Interface, err error) func() (gene.Interface, error) {
	return func() (gene.Interface, error) {
		tmpG, tmpE := g, err
		g, err = nil, io.EOF
		return tmpG, tmpE
	}
}

// Test ReadFromFunc
var readFromFuncTests = []struct {
	name string
	f    func() (gene.Interface, error)
	g    gene.Interface
	err  string
}{
	{
		f:   geneGenerator(&gene.Gene{ID: "A"}, nil),
		g:   &gene.Gene{ID: "A"},
		err: "",
	},
	{
		f:   geneGenerator(nil, errors.New("returning error")),
		g:   nil,
		err: "returning error",
	},
	{
		f:   geneGenerator(&gene.Gene{ID: "B"}, io.EOF),
		g:   &gene.Gene{ID: "B"},
		err: "",
	},
}

func TestReadFromFunc(t *testing.T) {
	for _, tt := range readFromFuncTests {
		sc := NewScannerFromFunc(tt.f)

		ok := sc.Next()
		g, err := sc.Gene(), sc.Error()

		// Next() should not return true when there is an error.
		if ok && err != nil {
			t.Errorf("%s: unexpected error %v", tt.name, err)
		}

		if tt.err != "" {
			if err == nil || !strings.Contains(err.Error(), tt.err) {
				t.Errorf("%s: error %q, want error %q", tt.name, err, tt.err)
			}
			continue
		} else if err != nil {
			t.Errorf("%s: unexpected error %v", tt.name, err)
			continue
		}

		if !reflect.DeepEqual(g, tt.g) {
			t.Errorf("%s: out gene=%v want %v", tt.name, g, tt.g)
		}

		// Second call to Next() should return false.
		if sc.Next() {
			t.Errorf("last cast call to next returned true instead of false")
		}
	}
}

// Test GeneReader's Read
var readTests = []struct {
	Name             string
	Input            string
	Sorted           bool
	Error            string
	GeneCnt, FeatCnt int
	IDs, Chrs        []string
	Starts, Ends     []int
	Orientations     []feat.Orientation
}{
	{
		Name: "Mixed sorted",
		Input: "" +
			"X\t.\texon\t10\t20\t0\t+\t.\tgene_id A; transcript_id A1;\n" +
			"X\t.\texon\t15\t30\t0\t+\t.\tgene_id A; transcript_id A2;\n" +
			"X\t.\texon\t20\t45\t0\t+\t.\tgene_id A; transcript_id A3;\n" +
			"X\t.\texon\t25\t30\t0\t+\t.\tgene_id A; transcript_id A4;\n" +
			"Y\t.\texon\t50\t90\t0\t-\t.\tgene_id B; transcript_id B1;\n" +
			"Y\t.\tfoo\t100\t200\t0\t-\t.\tgene_id B; transcript_id B1;\n" +
			"Y\t.\tstop_codon\t60\t62\t0\t-\t.\tgene_id B; transcript_id B1;\n" +
			"Y\t.\tstart_codon\t71\t73\t0\t-\t.\tgene_id B; transcript_id B1;\n",
		Sorted:       true,
		GeneCnt:      2,
		FeatCnt:      5,
		IDs:          []string{"A", "B"},
		Chrs:         []string{"X", "Y"},
		Starts:       []int{9, 49},
		Ends:         []int{45, 90},
		Orientations: []feat.Orientation{feat.Forward, feat.Reverse},
	},
	{
		Name: "Single forward coding gene",
		Input: "" +
			"X\t.\texon\t2\t70\t0\t+\t.\tgene_id C; transcript_id C1;\n" +
			"X\t.\texon\t80\t90\t0\t+\t.\tgene_id C; transcript_id C1;\n" +
			"X\t.\tstart_codon\t60\t62\t0\t+\t.\tgene_id C; transcript_id C1;\n" +
			"X\t.\tstop_codon\t81\t83\t0\t+\t.\tgene_id C; transcript_id C1;\n",
		Sorted:       true,
		GeneCnt:      1,
		FeatCnt:      1,
		IDs:          []string{"C"},
		Chrs:         []string{"X"},
		Starts:       []int{1},
		Ends:         []int{90},
		Orientations: []feat.Orientation{feat.Forward},
	},
	{
		Name: "Single reverse coding gene",
		Input: "" +
			"X\t.\texon\t30\t50\t0\t-\t.\tgene_id D; transcript_id D1;\n" +
			"X\t.\texon\t80\t99\t0\t-\t.\tgene_id D; transcript_id D1;\n" +
			"X\t.\tstop_codon\t40\t42\t0\t-\t.\tgene_id D; transcript_id D1;\n" +
			"X\t.\tstart_codon\t91\t93\t0\t-\t.\tgene_id D; transcript_id D1;\n",
		Sorted:       true,
		GeneCnt:      1,
		FeatCnt:      1,
		IDs:          []string{"D"},
		Chrs:         []string{"X"},
		Starts:       []int{29},
		Ends:         []int{99},
		Orientations: []feat.Orientation{feat.Reverse},
	},
	{
		Name:   "Missing gene_id",
		Input:  "X\t.\texon\t30\t50\t0\t-\t.\ttranscript_id D1;\n",
		Sorted: true,
		Error:  "gff: empty grouping gene_id field",
	},
	{
		Name:   "Missing transcript_id",
		Input:  "X\t.\texon\t30\t50\t0\t-\t.\tgene_id D;\n",
		Sorted: true,
		Error:  "gff: empty grouping transcript_id field",
	},
	{
		Name: "Features on different location",
		Input: "" +
			"X\t.\texon\t30\t50\t0\t-\t.\tgene_id F; transcript_id F1;\n" +
			"Y\t.\texon\t80\t99\t0\t-\t.\tgene_id F; transcript_id F1;\n",
		Sorted: true,
		Error:  "geneio: features on varying location for gene F",
	},
	{
		Name: "Features on different orientation",
		Input: "" +
			"X\t.\texon\t30\t50\t0\t-\t.\tgene_id G; transcript_id G1;\n" +
			"X\t.\texon\t80\t99\t0\t+\t.\tgene_id G; transcript_id G1;\n",
		Sorted: true,
		Error:  "geneio: features with varying orientation for gene G",
	},
	{
		Name: "Single non coding gene with adjacent exons",
		Input: "" +
			"X\t.\texon\t30\t50\t0\t-\t.\tgene_id H; transcript_id H1;\n" +
			"X\t.\texon\t51\t99\t0\t-\t.\tgene_id H; transcript_id H1;\n",
		Sorted:       true,
		GeneCnt:      1,
		FeatCnt:      1,
		IDs:          []string{"H"},
		Chrs:         []string{"X"},
		Starts:       []int{29},
		Ends:         []int{99},
		Orientations: []feat.Orientation{feat.Reverse},
	},
	{
		Name: "Single non coding gene with overlapping exons",
		Input: "" +
			"X\t.\texon\t2\t70\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\texon\t70\t90\t0\t+\t.\tgene_id F; transcript_id F1;\n",
		Sorted: true,
		Error:  "exons overlap",
	},
	{
		Name: "Single coding gene with adjacent exons",
		Input: "" +
			"X\t.\texon\t10\t100\t0\t-\t.\tgene_id I; transcript_id I1;\n" +
			"X\t.\tstop_codon\t40\t42\t0\t-\t.\tgene_id I; transcript_id I1;\n" +
			"X\t.\tstart_codon\t91\t93\t0\t-\t.\tgene_id I; transcript_id I1;\n",
		Sorted:       true,
		GeneCnt:      1,
		FeatCnt:      1,
		IDs:          []string{"I"},
		Chrs:         []string{"X"},
		Starts:       []int{9},
		Ends:         []int{100},
		Orientations: []feat.Orientation{feat.Reverse},
	},
	{
		Name: "Single coding gene with overlapping exons",
		Input: "" +
			"X\t.\texon\t2\t70\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\texon\t70\t90\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\tstart_codon\t60\t62\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\tstop_codon\t81\t83\t0\t+\t.\tgene_id F; transcript_id F1;\n",
		Sorted: true,
		Error:  "exons overlap",
	},
	{
		Name: "Single coding gene with no stop_codon",
		Input: "" +
			"X\t.\texon\t2\t70\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\texon\t71\t90\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\tstart_codon\t60\t62\t0\t+\t.\tgene_id F; transcript_id F1;\n",
		Sorted: true,
		Error:  "geneio: only one of start/stop codon found for F1",
	},
	{
		Name: "Single coding gene with no start_codon",
		Input: "" +
			"X\t.\texon\t2\t70\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\texon\t71\t90\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\tstop_codon\t81\t83\t0\t+\t.\tgene_id F; transcript_id F1;\n",
		Sorted: true,
		Error:  "geneio: only one of start/stop codon found for F1",
	},
}

func TestGeneReaderRead(t *testing.T) {
	for _, tt := range readTests {
		r := GeneReader{
			r: &FeatureReaderImpl{
				r:   gff.NewReader(strings.NewReader(tt.Input)),
				GID: "gene_id",
				TID: "transcript_id",
			},
		}

		genes, err := r.ReadAll()
		if tt.Error != "" {
			if err == nil || !strings.Contains(err.Error(), tt.Error) {
				t.Errorf("%s: error %q, want error %q", tt.Name, err, tt.Error)
			}
			continue
		} else if err != nil {
			t.Errorf("%s: unexpected error %v", tt.Name, err)
			continue
		}
		genesSorted := Genes(genes)
		sort.Sort(genesSorted)

		var feats []feat.Feature
		var ids, chrs []string
		var starts, ends []int
		var orientations []feat.Orientation
		for _, g := range genesSorted {
			feats = append(feats, g.Features()...)
			ids = append(ids, g.Name())
			starts = append(starts, g.Start())
			ends = append(ends, g.End())
			chrs = append(chrs, g.Location().Name())
			orientations = append(orientations, g.Orientation())
		}
		if len(genes) != tt.GeneCnt {
			t.Errorf("%s: out gene count=%d want %d", tt.Name, len(genes), tt.GeneCnt)
		}
		if len(feats) != tt.FeatCnt {
			t.Errorf("%s: out feat count=%d want %d", tt.Name, len(feats), tt.FeatCnt)
		}
		if !reflect.DeepEqual(ids, tt.IDs) {
			t.Errorf("%s: out ids=%q want %q", tt.Name, ids, tt.IDs)
		}
		if !reflect.DeepEqual(starts, tt.Starts) {
			t.Errorf("%s: out starts=%d want %d", tt.Name, starts, tt.Starts)
		}
		if !reflect.DeepEqual(ends, tt.Ends) {
			t.Errorf("%s: out ends=%d want %d", tt.Name, ends, tt.Ends)
		}
		if !reflect.DeepEqual(chrs, tt.Chrs) {
			t.Errorf("%s: out chrs=%q want %q", tt.Name, chrs, tt.Chrs)
		}
		if !reflect.DeepEqual(orientations, tt.Orientations) {
			t.Errorf("%s: out orientations=%q want %q", tt.Name, orientations, tt.Orientations)
		}
	}
}

func BenchmarkGeneReaderRead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		tt := readTests[0]
		r := GeneReader{
			r: &FeatureReaderImpl{
				r:   gff.NewReader(strings.NewReader(tt.Input)),
				GID: "gene_id",
				TID: "transcript_id",
			},
		}
		_, _ = r.ReadAll()
	}
}
