-- Verify model_model

BEGIN;

SELECT genome_model_id, id, name, sample_name, processing_profile_id,
    data_directory, comparable_normal_model_id, current_running_build_id,
    last_complete_build_id, user_name, creation_date, auto_assign_inst_data,
    auto_build_alignments, subject_id, subject_class_name, keep,
    build_granularity_unit, build_granularity_value, limit_inputs_to_id,
    is_default, subclass_name,auto_build, auto_build, build_requested
FROM model.model
WHERE FALSE;

ROLLBACK;
