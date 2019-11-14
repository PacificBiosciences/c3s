// Author: Armin Töpfer

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>

#include <pbcopper/cli2/CLI.h>
#include <pbcopper/data/Read.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BgzipFastaWriter.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaWriter.h>
#include <pbbam/FastqReader.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>

#include "third-party/spoa/spoa.hpp"

namespace PacBio {
class Abort : public std::exception
{};

namespace OptionNames {
// clang-format off
const CLI_v2::Option FastaName{
R"({
    "names" : ["name"],
    "description" : "Name of the fasta record.",
    "type" : "string"
})"};
const CLI_v2::PositionalArgument Input {
R"({
    "name" : "input.bam|xml|fastq|fq",
    "description" : "Source BAM, DATASET, or FASTQ file"
})"};

const CLI_v2::PositionalArgument Output {
R"({
    "name" : "output.fa|fa.gz",
    "description" : "Output FASTA or FASTA BGZIP"
})"};

// clang-format on
}  // namespace OptionNames

struct Settings
{
    const std::string CLI;
    const std::vector<std::string> InputFiles;
    const std::string FastaName;

    // Parses the provided CLI_v2::Results and retrieves a defined set of options.
    Settings(const PacBio::CLI_v2::Results& options)
        : CLI(options.InputCommandLine())
        , InputFiles(options.PositionalArguments())
        , FastaName(options[OptionNames::FastaName])
    {}

    // Given the description of the tool and its version, create all
    // necessary CLI_v2::Options for the victor executable.
    static PacBio::CLI_v2::Interface CreateCLI()
    {
        PacBio::CLI_v2::Interface i{"c3s", "Generate consensus of CCS reads", "0.0.1"};
        i.DisableNumThreadsOption();

        // clang-format off
        i.AddPositionalArguments({
            OptionNames::Input,
            OptionNames::Output
        });
        i.AddOption({
            OptionNames::FastaName
        });
        // clang-format on

        return i;
    }
};

struct SimpleBamParser
{
    static std::unique_ptr<BAM::internal::IQuery> BamQuery(const std::string& filePath)
    {
        BAM::DataSet ds(filePath);
        const auto filter = BAM::PbiFilter::FromDataSet(ds);
        std::unique_ptr<BAM::internal::IQuery> query(nullptr);
        if (filter.IsEmpty())
            query.reset(new BAM::EntireFileQuery(ds));
        else
            query.reset(new BAM::PbiFilterQuery(filter, ds));
        return query;
    }
};

struct SimpleRead
{
    std::string Seq;
    std::string Qual;
};

struct Workflow
{
    static int Runner(const PacBio::CLI_v2::Results& options)
    {
        Settings settings{options};
        if (settings.InputFiles.size() != 2) {
            PBLOG_FATAL << "Please provide the unaligned BAM file and output BAM file!";
            throw Abort();
        }
        const auto input = settings.InputFiles[0];
        if (!PacBio::Utility::FileExists(input)) {
            PBLOG_FATAL << "Input data file does not exist: " << input;
            throw Abort();
        }

        const auto output = settings.InputFiles[1];
        if (PacBio::Utility::FileExists(output))
            PBLOG_WARN << "Overwriting existing output file: " << output;

        const std::vector<std::string> outputFormats{".fa", ".fa.gz", ".fasta", ".fasta.gz"};
        if (std::all_of(
                outputFormats.cbegin(), outputFormats.cend(),
                [&output](const std::string& end) { return !boost::iends_with(output, end); }))
            PBLOG_WARN << "Unknown output file type. Please refer to --help!";

        std::vector<SimpleRead> reads;
        if (boost::iends_with(input, ".bam") || boost::iends_with(input, ".xml")) {
            auto query = SimpleBamParser::BamQuery(input);
            for (auto& read : *query) {
                if (read.IsMapped()) {
                    PBLOG_FATAL << "Input records must be unaligned. Offending record: "
                                << read.FullName();
                    throw Abort();
                }
                reads.push_back({read.Sequence(), read.Qualities().Fastq()});
            }
        } else if (boost::iends_with(input, ".fq") || boost::iends_with(input, ".fastq") ||
                   boost::iends_with(input, ".fq.gz") || boost::iends_with(input, ".fastq.gz")) {
            for (auto& read : BAM::FastqReader{input}) {
                const std::string seq = read.Bases();
                const std::string qual = read.Qualities().Fastq();
                if (seq.size() != qual.size()) {
                    PBLOG_FATAL << "Sequence and qualities must have same length! Offending read: "
                                << read.Name();
                    throw Abort();
                }
                reads.push_back({std::move(seq), std::move(qual)});
            }
        } else {
            PBLOG_FATAL << "Unknown input file type. Please refer to --help!";
            throw Abort();
        }

        auto alignment_engine = spoa::createAlignmentEngine(spoa::AlignmentType::kSW, 1, -2, -2);
        auto graph = spoa::createGraph();
        graph->add_alignment(spoa::Alignment(), reads[0].Seq, reads[0].Qual);

        const int32_t numReads = reads.size();
        for (int32_t i = 1; i < numReads; ++i) {
            spoa::Alignment alignment =
                alignment_engine->align_sequence_with_graph(reads[i].Seq, graph);
            graph->add_alignment(alignment, reads[i].Seq, reads[i].Qual);
        }

        const std::string consensus = graph->generate_consensus(std::max(numReads / 2, 1));

        std::unique_ptr<BAM::IFastaWriter> writer;
        if (boost::iends_with(output, ".gz"))
            writer = std::make_unique<BAM::BgzipFastaWriter>(output);
        else
            writer = std::make_unique<BAM::FastaWriter>(output);
        std::string fastaName{settings.FastaName};
        if (fastaName.empty()) fastaName = "consensus";
        writer->Write(fastaName, consensus);

        return EXIT_SUCCESS;
    }
};
}  // namespace PacBio

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, PacBio::Settings::CreateCLI(),
                                   &PacBio::Workflow::Runner);
    } catch (const PacBio::Abort& e) {
        return EXIT_FAILURE;
    } catch (const std::runtime_error& e) {
        std::cerr << "ERROR: " << e.what();
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
