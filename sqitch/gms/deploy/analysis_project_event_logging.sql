-- Deploy analysis_project_event_logging
-- requires: timeline_base_permission

BEGIN;

  CREATE TABLE IF NOT EXISTS timeline.analysis_project_event_type (
    id character varying NOT NULL,
    CONSTRAINT anp_event_type_pk PRIMARY KEY(id)
  );

  REVOKE ALL ON TABLE timeline.analysis_project_event_type FROM PUBLIC;
  GRANT ALL ON TABLE timeline.analysis_project_event_type TO genome;
  GRANT SELECT,INSERT,UPDATE,DELETE ON TABLE timeline.analysis_project_event_type TO "gms-user";

  CREATE TABLE IF NOT EXISTS timeline.analysis_project (
      status character varying(255) NOT NULL,
      is_cle boolean NOT NULL,
      run_as character varying(64) NOT NULL,
      CONSTRAINT anp_event_anp_fk FOREIGN KEY (object_id) REFERENCES config.analysis_project(id) MATCH FULL,
      CONSTRAINT anp_event_type_fk FOREIGN KEY (name) REFERENCES timeline.analysis_project_event_type(id) MATCH FULL
  ) INHERITS (timeline.base);

  CREATE INDEX timeline_anp_object_id_index ON timeline.analysis_project (object_id);

  REVOKE ALL ON TABLE timeline.analysis_project FROM PUBLIC;
  GRANT ALL ON TABLE timeline.analysis_project TO genome;
  GRANT SELECT,INSERT,UPDATE,DELETE ON TABLE timeline.analysis_project TO "gms-user";

COMMIT;
