/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.23
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "efg-locate"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "efg-locate"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0.0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *full_help_help; /**< @brief Print help, including hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * ignore_chars_arg;	/**< @brief Ignore these characters for the indexability property/pattern matching, breaking up each pattern into maximal strings of non-ignore characters.  */
  char * ignore_chars_orig;	/**< @brief Ignore these characters for the indexability property/pattern matching, breaking up each pattern into maximal strings of non-ignore characters original value given at command line.  */
  const char *ignore_chars_help; /**< @brief Ignore these characters for the indexability property/pattern matching, breaking up each pattern into maximal strings of non-ignore characters help description.  */
  int approximate_flag;	/**< @brief Approximate pattern matching by greedily matching the pattern in the graph and starting over when the matching fails; output only the recognized matches spanning at least a full node (default=off).  */
  const char *approximate_help; /**< @brief Approximate pattern matching by greedily matching the pattern in the graph and starting over when the matching fails; output only the recognized matches spanning at least a full node help description.  */
  int approximate_edge_match_min_count_arg;	/**< @brief Consider any approximate occurrence valid if the pattern substring occurs at most COUNT times in the edges (default='0').  */
  char * approximate_edge_match_min_count_orig;	/**< @brief Consider any approximate occurrence valid if the pattern substring occurs at most COUNT times in the edges original value given at command line.  */
  const char *approximate_edge_match_min_count_help; /**< @brief Consider any approximate occurrence valid if the pattern substring occurs at most COUNT times in the edges help description.  */
  int approximate_edge_match_min_count_heuristic_flag;	/**< @brief Together with the --approximate-edge-match-min-count option, compute the edge occurrences occurring at most COUNT times only after not finding any long match of some read substring in the graph spanning 3+ nodes (default=off).  */
  const char *approximate_edge_match_min_count_heuristic_help; /**< @brief Together with the --approximate-edge-match-min-count option, compute the edge occurrences occurring at most COUNT times only after not finding any long match of some read substring in the graph spanning 3+ nodes help description.  */
  int approximate_min_coverage_arg;	/**< @brief Consider approximate occurrences as valid if they cover at least PERC % of the pattern (default='0').  */
  char * approximate_min_coverage_orig;	/**< @brief Consider approximate occurrences as valid if they cover at least PERC % of the pattern original value given at command line.  */
  const char *approximate_min_coverage_help; /**< @brief Consider approximate occurrences as valid if they cover at least PERC % of the pattern help description.  */
  int approximate_stats_flag;	/**< @brief Output statistics for each read in stdout (default=off).  */
  const char *approximate_stats_help; /**< @brief Output statistics for each read in stdout help description.  */
  int reverse_complement_flag;	/**< @brief Match also the reverse complement of the patterns and output the results as a reverse graph path (default=off).  */
  const char *reverse_complement_help; /**< @brief Match also the reverse complement of the patterns and output the results as a reverse graph path help description.  */
  int rename_reverse_complement_flag;	/**< @brief When matching the reverse complement of patterns, consider them as a distinct patterns by prepending 'rev_' to its name (default=off).  */
  const char *rename_reverse_complement_help; /**< @brief When matching the reverse complement of patterns, consider them as a distinct patterns by prepending 'rev_' to its name help description.  */
  int split_output_matches_flag;	/**< @brief In approximate mode (--approximate), split long matches into node matches (default=off).  */
  const char *split_output_matches_help; /**< @brief In approximate mode (--approximate), split long matches into node matches help description.  */
  int split_output_matches_graphaligner_flag;	/**< @brief Same as --split-output-matches, but filter out node matches of length 1 (for use with GraphAligner --extend) (default=off).  */
  const char *split_output_matches_graphaligner_help; /**< @brief Same as --split-output-matches, but filter out node matches of length 1 (for use with GraphAligner --extend) help description.  */
  long threads_arg;	/**< @brief Number of compute threads (default='-1').  */
  char * threads_orig;	/**< @brief Number of compute threads original value given at command line.  */
  const char *threads_help; /**< @brief Number of compute threads help description.  */
  int overwrite_flag;	/**< @brief Overwrite the output file, if it exists (default=off).  */
  const char *overwrite_help; /**< @brief Overwrite the output file, if it exists help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int full_help_given ;	/**< @brief Whether full-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int ignore_chars_given ;	/**< @brief Whether ignore-chars was given.  */
  unsigned int approximate_given ;	/**< @brief Whether approximate was given.  */
  unsigned int approximate_edge_match_min_count_given ;	/**< @brief Whether approximate-edge-match-min-count was given.  */
  unsigned int approximate_edge_match_min_count_heuristic_given ;	/**< @brief Whether approximate-edge-match-min-count-heuristic was given.  */
  unsigned int approximate_min_coverage_given ;	/**< @brief Whether approximate-min-coverage was given.  */
  unsigned int approximate_stats_given ;	/**< @brief Whether approximate-stats was given.  */
  unsigned int reverse_complement_given ;	/**< @brief Whether reverse-complement was given.  */
  unsigned int rename_reverse_complement_given ;	/**< @brief Whether rename-reverse-complement was given.  */
  unsigned int split_output_matches_given ;	/**< @brief Whether split-output-matches was given.  */
  unsigned int split_output_matches_graphaligner_given ;	/**< @brief Whether split-output-matches-graphaligner was given.  */
  unsigned int threads_given ;	/**< @brief Whether threads was given.  */
  unsigned int overwrite_given ;	/**< @brief Whether overwrite was given.  */

  char **inputs ; /**< @brief unnamed options (options without names) */
  unsigned inputs_num ; /**< @brief unnamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];
/** @brief all the lines making the full help output (including hidden options) */
extern const char *gengetopt_args_info_full_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the full help (including hidden options)
 */
void cmdline_parser_print_full_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
