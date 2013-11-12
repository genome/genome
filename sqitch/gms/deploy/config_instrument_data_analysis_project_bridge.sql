-- Deploy config_instrument_data_analysis_project_bridge
-- requires: instrument_data
-- requires: config_analysis_project

BEGIN;

CREATE TABLE config.instrument_data_analysis_project_bridge (
  id character varying(64) PRIMARY KEY,
  instrument_data_id character varying(64) NOT NULL REFERENCES instrument.data(id),
  analysis_project_id character varying(64) NOT NULL REFERENCES config.analysis_project(id),
  created_at timestamp(6) NOT NULL,
  updated_at timestamp(6) NOT NULL,
  status character varying NOT NULL,
  reason character varying,
  fail_count integer NOT NULL,
  UNIQUE(instrument_data_id, analysis_project_id)
);

COMMIT;
