-- Deploy timeline_analysis_project_event_type_values
-- requires: analysis_project_event_logging

BEGIN;

-- It's written this way because these validation values are required for
-- the system to function, and were already in the table, but there was no
-- tracking for how they should get into the table in the first place

DO $$
    DECLARE vals varchar[] := ARRAY['config_added','status_changed','model_created',
                                    'cle_changed','instrument_data_assigned'];
    DECLARE i varchar;
BEGIN
    FOREACH i IN ARRAY vals
    LOOP
        BEGIN
            insert into timeline.analysis_project_event_type values (i);
        EXCEPTION when unique_violation THEN
            -- do nothing
        END;
    END LOOP;
END;
$$ language 'plpgsql';

COMMIT;
